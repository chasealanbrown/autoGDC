

def metadata_json_to_df(json_filepath: str,
                        null_type_strings: list = ["N/A"])-> pd.DataFrame:
  """
  Summary:
    Converting the full GDC metadata from json format to a Dataframe

  Notes:
    The typical JSON structure returned by GDC
      (considering only nested data structures)

      data:
        -cases
          - samples
            - portions
              - analytes
                - aliquots
          - follow_ups
          - demographic
          - diagnosis
        -analysis
  """
  # First, convert the json to a dataframe
  with open(json_filepath) as f:
    # The field "data.hits" is of most interest
    df = pd.DataFrame(json.loads(f.read())["data"]["hits"])

  LOG.debug("Converting {json_filepath}\
            to dataframe - df prior to transforms:\n {df}")

  # Expand all of the dictionaries and lists within each column
  #   Neseted structure shown by indentation
  expand_order = [
          "cases",
            "samples",
              "portions",
                "analytes",
                  "aliquots",
            "follow_ups",
            "demographic",
            "diagnoses",
          "analysis",
          "archive", # ?
          "project", # ?
          ]
  # First set file UUID (`id`) to index, as the `column_expand()` will
  #   add more inidces as the lists and dictionaries in each cell are expanded
  df = df.set_index("id")
  for nested_col in expand_order:
    df = column_expand(df, nested_col)

  # Make sure index is the same every time
  #   this function is run on a JSON from GDC's REST API
  #   (fill missing indices with zeros)
  index_columns = [
          "id",
          "nested_cases_index",
            "nested_samples_index",
              "nested_portions_index",
                "nested_analytes_index",
                  "nested_aliquots_index",
            "nested_follow_ups_index",
            "nested_diagnoses_index",
          ]
  # Use the index of `df` to check what is missing
  missing_index_cols = set(index_columns).difference(set(df.index.names))
  # Now reset in order to add back in the missing columns
  df = df.reset_index()
  for missing_col in list(missing_index_cols):
    df[missing_col] = 0

  # Now we set the index again, with more certainty that it is correct
  #   (The file UUID and these nested indices are necessary
  #    to join with metadata database)
  df = df.set_index(index_columns)

  # Convert strings denoting nulls into proper null values
  re_match_fullstr = [f"^{str(s).lower()}$" for s in null_type_strings]
  replacements = {'|'.join(re_match_fullstr): None}
  df = df.replace(replacements, regex=True)

  # Drop columns containing NaN data
  #   (Only those which have absolutely no data)
  df = df.dropna(axis = 1, how = "all")

  # Age is in days - change it to years
  if "age_at_diagnosis" in df.columns:
    df["age"] = np.round(df["age_at_diagnosis"]/365, 2)

  LOG.debug(f"Converted json file to dataframe with {df.shape[0]} rows.")
  return df


def column_expand(df, col: str):
  df = df.sort_index()
  try:
    colvals = df[col].dropna().tolist()

  except KeyError as e:
    LOG.debug(e)
    LOG.debug(f"DataFrame does not appear to have column: {col}")
    return df

  try:
    firstelem = colvals[0]
  except KeyError as e:
    LOG.debug(e)
    LOG.debug("First element cannot be accessed")
    return df

  if df.index.names == [None]:
    df.index.names = ["index"]
  try:
    if isinstance(firstelem, list):

      # This expands columns with lists as elements
      expand_df = df.explode(col).sort_index()

      # Since `.explode()` duplicates the index, and doesn't provide
      #   a multi-index, we create one here
      # Nest count is essentially how many items with in each list of each elem
      nest_cnt = expand_df.groupby(level = expand_df.index.names)[col].count().sort_index()
      assert len(nest_cnt) == len(df)

      # Make sure that any NaN values are still counted such that the length of
      #   the dateframe will end up being the same size
      nest_cnt[nest_cnt == 0] = 1
      assert nest_cnt.sum() == len(expand_df)

      # We then create a new index, where the nested count is rolled out
      #   and concatenated into one list
      #   (i.e. a multiindex with a counter for duplicates)
      muix_vals = np.array(list(iterchain(*map(range, nest_cnt.values))))
      expand_df[f"nested_{col}_index"] = muix_vals
      expand_df.set_index(f"nested_{col}_index", append = True, inplace = True)

      # Now recursively look for lists within this column
      return column_expand(expand_df, col)

    if isinstance(firstelem, dict):
      nan_replaced_series = df[col].where(df[col].notna(), lambda x: [{}])
      tmp_df = pd.DataFrame(nan_replaced_series.values.tolist(),
                            index = df.index)

    else:
      LOG.debug(f"The column {col} type is not yet supported")
      return df

    # Rename columns if necessary
    to_rename = set(df.columns).intersection(set(tmp_df.columns))
    prefix = f"{col}_"
    tmp_df.rename(columns = {c:f"{prefix}{c}" for c in tmp_df.columns if c in to_rename}, inplace = True)

    # Create full dataframe after expansion
    df = pd.concat([df, tmp_df], axis = 1).drop(col, axis =1).set_index(df.index)

  except Exception as e:
    LOG.warn(f"WARNING - There was an error in processing the column {col}")
    LOG.warn(e)

  return df


#def read_and_filter(filepath: str,
#                    subset_features: list = None):
#  """
#  Summary:
#    Reads a tsv file, formats any gene names, and renames it to GDC's `file_id`
#
#  Arguments:
#    filepath:
#      Path to tsv downloaded from GDC
#
#    subet_features:
#      Subset of the features as a list
#  """
#
#  # Read the series
#  #   (which was saved and compress during download phase)
#  series = pd.read_csv(filepath,
#                       sep = "\t",
#                       index_col = 0,
#                       header = None,
#                       engine = "c",
#                       squeeze = True)
#
#  # Rename it to be the filename / file_id
#  file_id = path.splitext(path.basename(filepath.rstrip(".gz")))[0]
#  series.name = file_id
#
#  series = series[~series.index.duplicated(keep="first")]
#
#  # Select subset of features
#  if subset_features:
#    series = series.reindex(subset_featuees)
#  return series


#@memoize
#def multifile_df(file_paths: list,
#                 subset_features: list = None):
#  """
#  Summary:
#    Reads a list of tsv files, formats any gene names,
#      and renames it to GDC's `file_id`
#
#  Arguments:
#    filepaths:
#      List of paths to tsv downloaded from GDC
#
#    subet_features:
#      Subset of the features as a list
#  """
#
##  kwargs = {"subset_features": subset_features}
#
#  # IO bound
#  # Parallel imap with progress bar
#  # series_iter = p_imap(partial(read_and_filter, **kwargs), file_paths)
#  series_list = [read_and_filter(fp, subset_features)
#                      for fp in tqdm(file_paths)]
#
#  # Concatenate all series
#  df = pd.DataFrame({s.name:s for s in series_list})
#
#  df.columns.name = "file_id"
#  return df.dropna(how = "all")


def contains_all_substrings(main_list, substr_list):
  """
  summary:
    determine if all substrings within a list of substrings are represented
      within a large list of strings

  arguments:
    main_list:
      list of strings that will be tested to ensure that all substrings are
        represented within this list.

    substr_list:
      list of substrings that need to be represented to return `true`

  example:
    substr_list = ["rna", "dna_methylation"]
    main_list = ["file_1_rna", "file_1_dna_methylation", "file_1_rna_2"]
    constains_all_substrings(main_list, substr_list)
    # >>> true

    main_list = ["file_1_rna", "file_1_rna_2"]
    contains_all_substrings(main_list, substr_list)
    # >>> false
  """
  return all(True if any((substr in fullstr for fullstr in main_list))
             else False
             for substr in substr_list)


def subset_paired_assay(mdf):
  """
  Summary:
    Returns a subset of the metadata dataframe wherein all of the cases
      contain all of the biological assays considered for the study

  Arguments:
    mdf:
      Metadata dataframe from the autoGDC object
  """
  # Get all of the assay types:
  assay_types = mdf["workflow_type"].unique()

  # Only select cases for which have paired assays
  mdf_file_lists = mdf.groupby("case_id")["workflow_type"].apply(list)
  mdf_paired_case_bool = mdf_file_lists.apply(
                            lambda x: contains_all_substrings(x, assay_types))
  paired_cases = mdf_paired_case_bool[mdf_paired_case_bool].index

  mdf = mdf[mdf["case_id"].isin(paired_cases)]
  return mdf


def quantile_normalize(df):
    """
    Summary:
      Performs quantile normalization on a dataframe of samples x features
        (rows x columns)
	Source:
      https://cmdlinetips.com/2020/06/computing-quantile-normalization-in-python/
    """
    # Input assumes samples as rows
    #   We transpose s.t. the implementation below
    #   (for features as rows) works
    df = df.T
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis = 0),
                             index=df.index,
                             columns=df.columns)
    df_mean = df_sorted.mean(axis = 1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn = df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    # Transpose back
    df_qn = df_qn.T
    return df_qn


#@memoize
def combined_region_collapsed_frame(dfs_dict,
                                   main_seq_data,
                                   seq_feature_metadata,
                                   region,
                                   agg_func = list,
                                   pos_col = "Position_to_TSS",
                                   max_seq_len = 20,
                                   value_name = "data_value",
                                   collapsed_mix_name = "RNA",
                                   sequence_mix_name = "DNA Methylation"):
    """
    Example Usage:

        dfs_dict = {"word_count": paired_df, "sentiment_score": paired_df2}
        combined_region_collapsed_frame(dfs_dict,
                                        main_seq_data = docs,
                                        seq_feature_metadata=feature_meta_info,
                                        region = "region",
                                        agg_func = np.median,
                                        pos_col = "position_in_region",
                                        max_seq_len = 4,
                                        value_name = "data_value")

    """

    all_region_df = combine_different_measures(dfs_dict)

    seq_region_collapsed = seq_collapsing(main_seq_data = main_seq_data,
                               seq_feature_metadata = seq_feature_metadata,
                               region = region,
                               agg_func = agg_func,
                               pos_col = pos_col,
                               max_seq_len = max_seq_len,
                               value_name = value_name)

    # Make multiindex columns
    all_region_df.columns = pd.MultiIndex.from_tuples(
                                            ((np.nan, col)
                                             for col in all_region_df.columns))

    LOG.info("Combining DNAm and RNA expression into one dataframe...");
    t0 = time.time()

    main_df = pd.concat([all_region_df, seq_region_collapsed],
                        keys = [collapsed_mix_name, sequence_mix_name],
                        axis = 1,
                        join = "inner")

    duration = round(t0 - time.time(), 1)

    LOG.info(f"Took {duration} seconds to concatenate the DNAm & RNA frames")
    return main_df

#@memoize
def combine_different_measures(df_dict):
    """
    Assuming a list of dataframes in (feature_rows x sample_cols) format,
        stack and order the dataframe by sample,
        then combine several measurements from different dataframes

    Examples:
        1) Multiple expression measurements (counts, FKPM, FKPM-UQ, user_normalized)
               in several different dataframes for many samples,
               combined to one dataframe.
            I.e. -- "region" could be a gene, "documents" could be samples,
                 and "word_count" could be rna-counts or FKPM


        2) Multiple dataframes describing the word cound and sentiment score
               for multiple documents, within various regions of the text.


     input:
         A dictionary of dataframe name and dataframe:

         df_dict = {"sentiment_score":sentiment_score, "word_count": word_count}

         word_count:                                               sentiment_score:
                              document_0  document_1                                   document_0  document_1
             region                                                   region
             mid                       8          73                  mid                      61          60
             end                      73          89       +          end                      72           7
             start                    60          22                  start                    74          74
             intro                    60          50                  intro                    22          33
             title                    17          52                  title                    13          77


     output:
         out_df:
                                     word_count  sentiment_score
                        region
             document_0 end                  73               72
                        intro                60               22
                        mid                   8               61
                        start                60               74
                        title                17               13
             document_1 end                  89                7
                        intro                50               33
                        mid                  73               60
                        start                22               74
                        title                52               77

    Example Usage:

        dfs_dict = {"word_count": peared_df, "sentiment_score": peared_df2}
        combined_df = combine_different_measures(dfs_dict)

    """

    # Remove null dataframes
    df_dict = {k:v for k,v in df_dict.items() if len(v)>0}

    # Take mean value for duplicated samples
    for name, df in df_dict.items():
        df_dict[name] = df.groupby(df.columns, axis=1).agg(np.mean)

    LOG.info("Combining all RNA documents to one dataframe...")
    t0 = time.time()

    # Create combined dataframes by stacking samples
    #   while new columns are appended for different measured variables
    output = pd.concat([v for v in df_dict.values()],
                       keys = [name for name in df_dict.keys()],
                       axis = 1)\
                    .stack()\
                    .swaplevel(0, 1, axis = 0)\
                    .sort_index(level = 0, axis = 0)


    rand_key = next(iter(df_dict.keys()))
    original_ix_name = df_dict[rand_key].index.name
    original_col_name = df_dict[rand_key].columns.name


    if original_col_name is None:
        original_col_name = "document"

    output.index.set_names([original_col_name, original_ix_name], inplace = True)

    duration = round(time.time() - t0, 1)
    LOG.info(f"Took {duration} seconds to combine the measures")
    return output


#@memoize
def seq_collapsing(main_seq_data,
                   seq_feature_metadata,
                   region,
#                    other_features,
                   agg_func = list, # np.median
                   pos_col = "Position_to_TSS",
                   max_seq_len = 20,
                   value_name = "data_value",
                   duplicate_agg_func = np.mean):
    """
    Assuming a list of dataframes in (feature_rows x sample_cols) format,
        stack and order the dataframe by sample,
        then combine several measurements from different dataframes

    Examples:
        1) Multiple expression measurements (counts, FKPM, FKPM-UQ, user_normalized)
               in several different dataframes for many samples,
               combined to one dataframe.
            I.e. -- "region" could be a gene, "documents" could be samples,
                 and "word_count" could be rna-counts or FKPM


        2) Multiple dataframes describing the word cound and sentiment score
               for multiple documents, within various regions of the text.


     input:
         main_seq_data: pd.DataFrame
             A (feature_row x sample_col) dataframe.
             For example, this could be a (Illumina 450k array chip loci x samples) dataframe.

         seq_feature_metadata: pd.DataFrame
             A (feature_row x feature_metadata) dataframe.
             For example, this could be a (Illumina 450k array chip loci x metadata) dataframe,
                 wherein metadata is "Gene_Symbol" or "Position_to_TSS", describing the assoicated values for that locus


        main_seq_data:                                    seq_feature_metadata:
            main_seq_data  document_0  document_1             seq_feature_metadata region  position_in_region  other_uninteresting_info
            feat_index                                        feat_index
            0                    0.00        0.11             0                       end                5930                      3797
            1                    0.77        0.23             1                     title                6851                      2635
            2                    0.85        0.05     +       2                     intro                8660                      4577
            3                    0.36        0.52             3                     intro                4260                      3191
            4                    0.97        0.04             4                     title                6185                      4774
            5                    0.79        0.87             5                       mid                3105                      7156
            6                    0.67        0.31             6                     intro                9559                      1077
            7                    0.00        0.39             7                       mid                  93                      5858
            8                    0.66        0.36             8                     intro                6085                       859
            9                    0.17        0.55             9                     start                3687                      7554


        region:
            This is the column in the metadata to collapse on.
            i.e. "Gene_symbol"

        other_features:
            Other features to be placed under each position index.
            i.e. "position in region" or "regulatory element"

     output:
         out_df:
            full_region_sub_ix           0                                         1 ...                                         3
                               data_values feat_index position_in_region data_values ... feat_index position_in_region data_values
            document   region
            document_0 end            0.63        2.0             5091.0        0.66 ...        1.0             9699.0         NaN
                       mid            0.95        4.0              423.0        0.20 ...        NaN                NaN         NaN
                       start          0.83        0.0             2891.0         NaN ...        NaN                NaN         NaN
            document_1 end            0.27        2.0             5091.0        0.23 ...        1.0             9699.0         NaN
                       mid            0.86        4.0              423.0        0.45 ...        NaN                NaN         NaN
                       start          0.95        0.0             2891.0         NaN ...        NaN                NaN         NaN
    """

    # Take mean value for duplicated samples
    main_seq_data = main_seq_data.groupby(main_seq_data.columns, axis=1).agg(duplicate_agg_func)

    if main_seq_data.index.name is None:
        feature_index = "feature_index"
    else:
        feature_index = main_seq_data.index.name

    if main_seq_data.columns.name is None:
        doc_col = "document"
    else:
        doc_col = main_seq_data.columns.name

    # Join the feature annotations with the data
    #   and stack the samples in order to obtain a very long series
    #
    # This returns a series with this format:
    #
    #     region  feat_index  position_in_region
    #     start   0           2891                document_0    0.83
    #                                             document_1    0.95
    #     end     1           9699                document_0    0.94
    #                                             document_1    0.05
    #             2           5091                document_0    0.63
    #                                             document_1    0.27
    #     mid     3           8801                document_0    0.20
    #                                             document_1    0.45
    #             4           423                 document_0    0.95
    #                                             document_1    0.86
    #     end     5           6018                document_0    0.66
    #                                             document_1    0.23
    #     dtype: float64
    #
    LOG.info("Processing and collapsing sequence data:\n"+\
                "Joining annotations and data..."); t0 = time.time()

    s = main_seq_data.join(seq_feature_metadata[[region, pos_col]])\
                        .reset_index()\
                        .set_index([region, pos_col, feature_index])\
                        .stack()

    LOG.info("Took "+str(round(time.time() - t0, 1)) + " seconds.\n"+\
              "Ordering by sample and gene position..."); t0 = time.time()

    # We then order by document/sample and the region/gene
    #
    # returning a dataframe of the following format:
    #
    #     document    region  position_in_region  feat_index
    #     document_0  end     5091                2             0.63
    #                         6018                5             0.66
    #                         9699                1             0.94
    #                 mid     423                 4             0.95
    #                         8801                3             0.20
    #                 start   2891                0             0.83
    #     document_1  end     5091                2             0.27
    #                         6018                5             0.23
    #                         9699                1             0.05
    #                 mid     423                 4             0.86
    #                         8801                3             0.45
    #                 start   2891                0             0.95
    #     Name: data_values, dtype: float64
    s.index.set_names([*s.index.names[:-1], doc_col], inplace = True)
    s.name = value_name
    s = s.reorder_levels(order = [s.index.names[-1], *s.index.names[:-1]])
    s = s.sort_index(level = s.index.names[:-1], axis = 0)


    LOG.info("Took "+str(round(time.time() - t0, 1)) + " seconds.\n"+\
              "Collapsing..."); t0 = time.time()

    if agg_func is not list:
        s = s.reset_index().pivot_table(index = [doc_col, region],
                                        values = [pos_col, value_name],
                                        aggfunc = agg_func)[value_name]
        s = s.to_frame()
        s.columns = pd.MultiIndex.from_tuples([(np.nan, c) for c in s.columns])

    else:
        # Make a multi-index with padding of the sequence
        #
        # This will be of the following form:
        #
        #                                                             document region  full_region_sub_ix
        #                     document   region full_region_sub_ix
        #                     document_0 end    0                   document_0    end                   0
        #                                       1                   document_0    end                   1
        #                                       2                   document_0    end                   2
        #                                       3                   document_0    end                   3
        #                                mid    0                   document_0    mid                   0
        #                                                 ...
        #                     document_1 end    0                   document_1    end                   0
        #                                                 ...
        #                                       3                   document_1    mid                   3
        #                                start  0                   document_1  start                   0
        #                                       1                   document_1  start                   1
        #                                       2                   document_1  start                   2
        #                                       3                   document_1  start                   3
        mux = pd.MultiIndex.from_product([s.index.get_level_values(doc_col).unique(),
                                          s.index.get_level_values(region).unique(),
                                          np.arange(max_seq_len)],
                                          names = [doc_col, region, "position_index"])

        # Remove position and feature index from index
        # i.e. keep document and region within the index (the index which is desired in the output)
        #
        # This will provide the following format of output:
        #
        #                                position_in_region  feat_index  data_values
        #             document   region
        #             document_0 end                   5091           2         0.63
        #                        end                   6018           5         0.66
        #                        end                   9699           1         0.94
        #                        mid                    423           4         0.95
        #                        mid                   8801           3         0.20
        #                        start                 2891           0         0.83
        #             document_1 end                   5091           2         0.27
        #                        end                   6018           5         0.23
        #                        end                   9699           1         0.05
        #                        mid                    423           4         0.86
        #                        mid                   8801           3         0.45
        #                        start                 2891           0         0.95
        s = s.reset_index(level=[2, 3])

        # Append the gene iteration index to the index
        #
        # This will provide the following format:
        #
        #                                  position_in_region  feat_index  data_values
        #             document   region
        #             document_0 end    0                5091           2         0.63
        #                               1                6018           5         0.66
        #                               2                9699           1         0.94
        #                        mid    0                 423           4         0.95
        #                               1                8801           3         0.20
        #                        start  0                2891           0         0.83
        #             document_1 end    0                5091           2         0.27
        #                               1                6018           5         0.23
        #                               2                9699           1         0.05
        #                        mid    0                 423           4         0.86
        #                               1                8801           3         0.45
        #                        start  0                2891           0         0.95
        s = s.set_index(s.groupby(level = [doc_col, region]).cumcount(), append = True)

        # Finally, reindex to the multi-index containing all of the sequence indices,
        # Then unstack these and sort each feature underneath it"s sequence index
        #
        # This provides the output formatted data
        s = s.reindex(mux).unstack().swaplevel(1, 0, axis = 1).sort_index(axis = 1)

        duration = round(t0-time.time(), 1)
        LOG.info(f"Took {duration} seconds to perform seq-collapse procedure.")

    return s
