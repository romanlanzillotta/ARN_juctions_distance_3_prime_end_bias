import pandas as pd


# files paths
junctions_file = 'isq_pcb1_final.txt'
sorted_gtf_file = 'isq_pcb1_sorted.gtf'
annotations_file = 'gencode.v42.annotation.gtf'
output_distances_file = 'distance.txt'

# column indexes (zero-based)
SORTED_GTF_ATTRIBUTE_GENE_ID_COLUMN = 0
SORTED_GTF_ATTRIBUTE_TRANSCRIPT_ID_COLUMN = 3
ANNOTATIONS_ATTRIBUTE_GENE_ID_COLUMN = 0
ANNOTATIONS_ATTRIBUTE_TRANSCRIPT_ID_COLUMN = 1

# junction_max_number. 0 = process all junctions in junctions file.
# use an integer > 0 if you want to limit the number of junctions processed.
junction_max_number = 0

asm_col_gene = SORTED_GTF_ATTRIBUTE_GENE_ID_COLUMN
asm_col_tr = SORTED_GTF_ATTRIBUTE_TRANSCRIPT_ID_COLUMN
ann_col_gene = ANNOTATIONS_ATTRIBUTE_GENE_ID_COLUMN
ann_col_tr = ANNOTATIONS_ATTRIBUTE_TRANSCRIPT_ID_COLUMN


# FUNCTIONS DEFINITIONS
def load_junctions(junctions_path):
    with open(junctions_path, 'r') as file:
        content = file.read()
        dict_ = eval(content)
        return dict_


def load_gtf(gtf_path):
    gtf_df = pd.read_csv(gtf_path, sep='\t',
                         comment='#', header=None,
                         names=['seqname', 'source', 'feature',
                                'start', 'end', 'score', 'strand',
                                'frame', 'attribute'])
    return gtf_df


def find_junction_rows(tuple_, gtf_df):
    """
    Returns the paired exon rows in the sorted GTF assembled file that form the junction.
    There may be more than one pair of rows.
    """
    filtered_df = gtf_df[(gtf_df['end'] == tuple_[0])
                         | (gtf_df['start'] == tuple_[1])]
    junction_indexes = []
    for i, chr_ in enumerate(filtered_df['seqname'].drop_duplicates()):
        df_chr = filtered_df[filtered_df['seqname'] == chr_]
        j = 0
        while j < (len(df_chr.index) - 1):
            if (df_chr.loc[df_chr.index[j], 'end'] == tuple_[0]) & (
                    df_chr.loc[df_chr.index[j + 1], 'start'] == tuple_[1]):
                junction_indexes.append(df_chr.index[j])
                junction_indexes.append(df_chr.index[j + 1])
            j += 1
    return filtered_df.loc[junction_indexes, :]


def fetch_gene_tr_id(tuple_, gtf_df):
    """
    Given start and end coordinates as input,
    fetches the gene id and transcript id from a sorted assembled gtf file.
    Returns a tuple.
    """
    gene_id = []
    gene_dict = {}
    filtered_df = find_junction_rows(tuple_, gtf_df)
    if filtered_df.shape[0] > 2:
        print(filtered_df)
    for _, row in filtered_df.iterrows():
        aux_gene = row['attribute'].split(";")[asm_col_gene]
        aux_gene = aux_gene.replace(r"gene_id", "").strip().strip(r'"')
        aux_transcript = row['attribute'].split(";")[asm_col_tr]
        aux_transcript = aux_transcript.replace(r"reference_transcript_id", "").strip().strip(r'"')
        gene_id.append(aux_gene)
        if not (aux_gene in gene_dict):
            gene_dict[aux_gene] = set()
        gene_dict[aux_gene].add(aux_transcript)

    if len(gene_dict) > 1:
        print("Warning: more than one gene ID found for junction (" +
              str(tuple_[0]) + ", " + str(tuple_[1]) + ").")
    return gene_id[0], gene_dict[gene_id[0]]


def fetch_genetable(gene_id, annog):
    """
    Fetches the table with all the features (genes, transcripts, exons, etc.)
    corresponding to the gene_id.
    Returns a pandas DataFrame with all the gene structure features.
    """
    gene = annog[annog['attribute'].str.contains(gene_id)]
    if len(gene) > 0:
        gene_start, gene_end = int(gene['start'].iloc[0]), int(gene['end'].iloc[0])
        seqname = gene['seqname'].iloc[0]
        # search for transcripts in the same chromosome within the gene range coords
        genetable = annotations[(annotations['seqname'] == seqname)
                                & (annotations['feature'].isin(['gene', 'transcript', 'exon']))
                                & (annotations['start'] >= gene_start)
                                & (annotations['end'] <= gene_end)]
        genetable = genetable[genetable['attribute'].str.contains(gene_id)]  # extra precaution
        gene_row = genetable[genetable['feature'] == 'gene']
        genetable = genetable[genetable['feature'] != 'gene']
        return gene_row, genetable
    else:
        return "", ""


def sum_exons(gene_table, idx_ini, idx_end):
    df = gene_table.loc[idx_ini:idx_end, 'start':'end']
    df.insert(df.shape[1], 'exon_length', df['end'] - df['start'] + 1)
    transcript_length = sum(df['exon_length'])
    return transcript_length


def create_dict_exons(gene_table):
    """
    Returns a dictionary with the range of rows that hold the exons for the given transcript index (key).

    transcript_idx: (exon_1_idx, exon_n_idx)
    """
    dict_ = {}
    idxs = list(gene_table.index)
    idxs_tr = list(gene_table[gene_table['feature'] == 'transcript'].index)
    for i in range(len(idxs_tr)):
        exon_1_idx = idxs[idxs.index(idxs_tr[i]) + 1]
        if i < (len(idxs_tr) - 1):
            exon_n_idx = idxs[idxs.index(idxs_tr[i + 1]) - 1]
        else:
            # last element
            exon_n_idx = idxs[-1]
        dict_[idxs_tr[i]] = (exon_1_idx, exon_n_idx)
    return dict_


def add_length_tr_col(tr_table, dict_exons, gene_table):
    """
    Return the same transcripts table with the calculated transcript length
    (sum of all of its exons) in the new column 'transcript_length'
    """
    transcripts_idx = list(tr_table.index)
    transcript_length = []
    for idx in transcripts_idx:
        idx_exons = dict_exons[idx]
        transcript_length.append(sum_exons(gene_table, idx_exons[0], idx_exons[1]))
    tr_table.insert(tr_table.shape[1], 'transcript_length', transcript_length)
    return tr_table


def fetch_transcripts_table(gene_table):
    """
    Receives a pandas DataFrame with a single gene structure features.
    Returns a pandas DataFrame with the transcripts table.
    """
    tr_table = gene_table[gene_table['feature'] == 'transcript']
    return tr_table


def extract_transcript_id(transcript_table):
    """
    Extracts the transcripts ids from annotation table,
    and adds them to the table as a column before the attributes.
    """
    transcript_id = []
    for _, row in transcript_table.iterrows():
        transcript_id.append(row['attribute'].split(";")[ann_col_tr].replace(r"transcript_id", "").strip().strip(r'"'))
    tbl = transcript_table
    tbl.insert(5, 'transcript id', transcript_id)
    tbl = tbl.iloc[:, 0:6]
    return tbl


def filter_transcripts_by_junction(junction, tbl_tr, rev_exon_dict, gtable):
    """
    This function returns all the possible transcripts from the given table tbl_tr in which the junction is present.
    If the junction coincides with any exon limits, then the transcripts that use those exons are returned.
    If the junction is not present in any annotated transcripts from the table, it is a novel junction.
        - If the novel junction is within any transcripts limits, those transcripts are returned.
        - If the novel junctions falls outside all transcripts limits, the transcripts with the closest
        exon are returned.
    """
    tbl_ex_tr = gtable[(gtable['end'] == junction[0])
                       | (gtable['feature'] == 'transcript')]
    tbl_ex_tr_idxs = list(tbl_ex_tr.index)
    tbl_ex_idxs = list(tbl_ex_tr[tbl_ex_tr['feature'] == 'exon'].index)
    if len(tbl_ex_idxs) > 0:
        # if the junction exists, the transcripts can be filtered
        # to those which contain exons separated by this junction
        tr_idx = []
        for i in tbl_ex_idxs:
            tr_idx.append(tbl_ex_tr_idxs[tbl_ex_tr_idxs.index(i) - 1])
        tbl_filtr = tbl_ex_tr.loc[tr_idx, :]
    else:
        # NOVEL JUNCTIONS
        # if it is a novel junction, we search the whole transcript table
        # the transcripts might not be annotated, but a requirement could be that the
        # junction at least overlaps or falls within the total transcript coords.
        tbl_filtr = tbl_tr[(tbl_tr['start'] <= junction[0]) & (tbl_tr['end'] >= junction[1])]
        if len(tbl_filtr) == 0:
            # the junction (start) falls completely or partially outside the transcript limits.
            # the transcript with the closest exon should be considered. There could be more than
            # one transcript with that exon. Only the outer exons are evaluated
            ends_exons = gtable[gtable.index.isin(rev_exon_dict.keys())]
            ends_exons.insert(ends_exons.shape[1], 'exon_dist_abs', abs(junction[0] - ends_exons['start'] + 1))
            # the closest exon might be used in more than one isoform and the minimum distance also might be the same
            # for several exons starting at the same position but with different lengths. We will consider them all.
            # with the minimum distance found the table is filtered again.
            min_distance = ends_exons['exon_dist_abs'].min()  # minimum absolute distance j start - exon start
            filtered_ends_exons = ends_exons[ends_exons['exon_dist_abs'] == min_distance]
            closest_exons_ids = list(filtered_ends_exons.index)
            transcripts_list = []
            for exon_idx in closest_exons_ids:
                transcripts_list.append(rev_exon_dict[exon_idx])
            tbl_filtr = tbl_tr[tbl_tr.index.isin(transcripts_list)]
        else:
            # One of two options: 1) return all transcripts where the junction falls within.
            # 2) Evaluate each exon. more time-consuming.
            # for now, 1)

            pass
    return tbl_filtr


def calc_longest_transcript(tbl_tr, lookup_exons, gene_table):
    """
    Returns a tuple for the longest transcript of the gene (idx_tr_table, transcript_id)
    Note: idx_tr_table = idx_gene_table
    """
    tbl_tr = add_length_tr_col(tbl_tr, lookup_exons, gene_table)
    id_max_length = tbl_tr['transcript_length'].idxmax()
    longest_tr = tbl_tr.loc[id_max_length, 'attribute'].split(";")[ann_col_tr]
    longest_tr = longest_tr.replace(r"transcript_id", "").strip().strip(r'"')
    return id_max_length, longest_tr


def fetch_exons(gtable, dict_exons, idx_tr):
    """
    Returns a DataFrame containing all the exons of the transcript.
    """
    exon_1_idx, exon_n_idx = dict_exons[idx_tr]
    df = gtable.loc[exon_1_idx:exon_n_idx, :]
    return df


def reverse_exon_dict(exon_dict):
    """
    Given a dictionary idx_tr: (exon_1_idx, exon_n_idx)
    Returns a dictionary exon_x_idx: idx_tr
    It contains all the firsts and lasts exons of each transcript
    """
    exons_rev_dict = {}
    for idx_transcript in exon_dict.keys():
        strart_end_exons = exon_dict[idx_transcript]
        for exon_idx in strart_end_exons:
            exons_rev_dict[exon_idx] = idx_transcript
    return exons_rev_dict


def calc_dist_3_end(junction, tr_exons):
    # find all exons with start >= junction 3' end
    complete_exons = tr_exons[(tr_exons['end'] > junction[1])
                              & (tr_exons['start'] >= junction[1])]
    indexes = list(complete_exons.index)
    if len(indexes) > 0:
        distance_3end = sum_exons(tr_exons, indexes[0], indexes[-1])
    else:
        # no exons were found after the splice junction end. This means that the splice junction is located
        # after all the exons of the transcript. In this case the distance should be taken fron the last end of exon
        # and will result in a negative value.
        distance_3end = tr_exons['end'].max() - junction[1] + 1

    # incomplete cases: the junction is de novo and overlaps an existing exon
    # complete overlap: junction falls within exon
    junction_within_exon = tr_exons[(tr_exons['start'] <= junction[0]) & (tr_exons['end'] >= junction[1])]
    if len(junction_within_exon) > 0:
        # print("junction_within_exon")
        # print(junction_within_exon.iloc[:,0:5])
        distance_3end = distance_3end + (junction_within_exon['end'].iloc[0] - junction[1] + 1)
    else:
        # partial overlap: junctions ends are shifted to a side of the nominal exon.
        # overlap at the 5' end of exon
        overlap_left = tr_exons[(tr_exons['start'] > junction[0]) & (tr_exons['start'] < junction[1])]
        if len(overlap_left) > 0:
            # print("overlap_left")
            # print(overlap_left)
            distance_3end = distance_3end + (overlap_left['end'].iloc[0] - junction[1] + 1)
        else:
            # overlap at the 3' end of exon
            overlap_right = tr_exons[(tr_exons['end'] > junction[0]) & (tr_exons['end'] < junction[1])]
            if len(overlap_right) > 0:
                distance_3end = distance_3end - (overlap_right['end'].iloc[0] - junction[0] + 1)
                # print("overlap_right")
                # print(overlap_right)
    return distance_3end


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    djunctions = load_junctions(junctions_file)  # load the dictionary file with junctions
    junctions = list(djunctions.keys())  # extract tuples (start, end)
    sorted_gtf = load_gtf(sorted_gtf_file)  # load the file with sorted assembled gtf
    annotations = load_gtf(annotations_file)  # load the annotations file

    if (type(junction_max_number) is int) and (junction_max_number > 0):
        junctions_limit = min(junction_max_number, len(junctions))
    else:
        junctions_limit = len(junctions)
    print("Total junctions to process:", junctions_limit)
    junctions_limit = 9
    with open(output_distances_file, "w") as output:
        out_string = "{"
        output.write(out_string)
        annogene = annotations[annotations['feature'] == 'gene']  # keep a reduced df version of gene rows
        i = 8
        junctions_limit -= 1
        gene_id, transcript_id = fetch_gene_tr_id(junctions[i],
                                                  sorted_gtf)  # the gene id needs to be fetched from the GTF file
        while i <= junctions_limit:
            if gene_id:
                gene_id_prev = gene_id
                gene_row, gene_table = fetch_genetable(gene_id,
                                                       annogene)  # the structure table from the annotation file,
                # for the gene id
                # print("gene_table")
                # print(gene_table.iloc[:,0:5])
                if len(gene_row) != 0:
                    lookup_exons_idx = create_dict_exons(gene_table)
                    tr_table = fetch_transcripts_table(gene_table)
                    reverse_lookup_exons = reverse_exon_dict(lookup_exons_idx)
                    # print(reverse_lookup_exons)
                    if len(tr_table) != 0:
                        while (i <= junctions_limit) and (gene_id_prev == gene_id):
                            # transcripts are filtered: if the junction is present in one or
                            # more transcript, keep those transcripts.
                            # if the junction is de novo: we keep the transcripts
                            # (one or many sharing the same exon)
                            # with the closest exon to the junction.
                            tr_table_filtr = filter_transcripts_by_junction(
                                junctions[i], tr_table,
                                reverse_lookup_exons,
                                gene_table)
                            # from the transcripts filtered, we keep the longest that contains
                            # the junction or if the novo,
                            # the longer one with the closest exon
                            idx_longest_tr, longest_tr_id = calc_longest_transcript(tr_table_filtr,
                                                                                    lookup_exons_idx, gene_table)

                            tr_ids_aux = transcript_id.copy()
                            if 'novel' in tr_ids_aux:
                                tr_ids_aux.remove('novel')

                            print("junction: ", i, junctions[i], gene_id, transcript_id)
                            if not (longest_tr_id in tr_ids_aux):
                                print("Assigned transcript(s) by assembler: ", transcript_id,
                                      "Longest transcript that contains the junction: ", longest_tr_id)
                                print(add_length_tr_col(extract_transcript_id(tr_table), lookup_exons_idx, gene_table))

                            # print("idx longest transcript:", idx_longest_tr, "transcript length:",
                            #      tr_table_filtr.loc[idx_longest_tr,'transcript_length'],
                            #      "transcript id:", longest_tr_id)

                            idx_longest_tr = gene_table[
                                (gene_table['feature'] == 'transcript') & gene_table['attribute'].str.contains(
                                    list(transcript_id)[0])].index[0]
                            print("idx transcript used:")
                            print(idx_longest_tr)
                            ex_table = fetch_exons(gene_table,
                                                   lookup_exons_idx,
                                                   idx_longest_tr)
                            # print("exon table:")
                            # print(ex_table.iloc[:,0:5])
                            distance = calc_dist_3_end(junctions[i], ex_table)
                            print("distance")
                            print(distance)
                            # if(distance<0):
                            #    print("junction after transcripts ends")
                            #    print(gene_table)
                            # write output to file
                            out_string = str(junctions[i]) + ": " + str(distance)
                            if i < junctions_limit:
                                out_string = out_string + ", "
                            output.write(out_string)
                            i += 1
                            if i <= junctions_limit:
                                gene_id, transcript_id = fetch_gene_tr_id(junctions[i],
                                                                          sorted_gtf)
                    else:
                        print("junction ", junctions[i],
                              ": Warning! there are no transcripts annotated in annotations file for gene id:",
                              gene_id, ". Skipped.")
                        i += 1
                        if i <= junctions_limit:
                            gene_id, transcript_id = fetch_gene_tr_id(junctions[i],
                                                                      sorted_gtf)
                else:
                    if gene_id.find("novel") != -1:
                        print("junction ", junctions[i],
                              ": Warning! gene id:", gene_id,
                              "is flagged as novel gene. There are no annotations for the transcripts. Skipped.")
                    else:
                        print("junction ", junctions[i],
                              ": Warning! gene id:", gene_id,
                              "not found in annotation file. Check that the annotations file version",
                              " provided is the same used for the assembled transcripts in the sorted GTF. Skipped.")
                    i += 1
                    if i <= junctions_limit:
                        gene_id, transcript_id = fetch_gene_tr_id(junctions[i],
                                                                  sorted_gtf)
            else:
                print("junction ", junctions[i],
                      ": Warning! gene id not found in sorted GTF file. Skipped.")
                i += 1
                if i <= junctions_limit:
                    gene_id, transcript_id = fetch_gene_tr_id(junctions[i],
                                                              sorted_gtf)
        output.write("}")

    print("Process finished.")
