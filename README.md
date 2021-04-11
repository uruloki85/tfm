# TFM

* Download data from GSE.
* Extract GSE accession and platform:

```bash
egrep "Series_platform_id|Series_geo_accession" GSE*_series_matrix.txt | cut -f2 | tr '\n' '\t' > series_and_platforms.csv
```

* Extract series matrix table:

```bash
sed -n '/^!series_matrix_table_begin/,/^!series_matrix_table_end/{p;/^!series_matrix_table_end/q}' GSE21501_series_matrix.txt | sed -e '1d' -e '$d'
```

The first `sed` cuts the content between the 2 patterns, `!series_matrix_table_begin` and `!series_matrix_table_end`, and exits when the second pattern is found.
The second `sed` removes the first and last line.

To apply this to all the files in the folder:

```bash
for file in *.txt; do
    sed -n '/^!series_matrix_table_begin/,/^!series_matrix_table_end/{p;/^!series_matrix_table_end/q}' "$file" | sed -e '1d' -e '$d' > "${file}_data"
done 
```

Result is saved in a new file with the suffix `_data`.
