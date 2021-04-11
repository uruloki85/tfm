#!/usr/bin/env python3.9
# -*- coding: UTF-8 -*-

import os
import argparse
import logging
import csv
# import pandas as pd
import GEOparse as geo


def read_series_and_platforms(filename):
    series_and_platforms = {}
    with open(filename, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            series_and_platforms[row[0]] = row[1]
        print(series_and_platforms)

    return series_and_platforms


def read_platform_data(platform, platform_dir):
    # dir = '/home/sabela/OneDrive/Personal/UOC/2020-21_02/TFM/data/platforms'
    platform_data = {}
    for entry in os.scandir(platform_dir):
        if not entry.is_file() or not entry.name.startswith(platform):
            continue

        path = entry.path
        with open(path, 'r') as f:
            csv_reader = csv.DictReader(f, delimiter='\t')
            for row in csv_reader:
                if not row:
                    continue  # Skip emtpy lines

                gene_id = row['ID']

                if 'GB_ACC' in row:
                    gene_bank_accession = row.get('GB_ACC')
                elif 'GB_LIST' in row:
                    gene_bank_accession = row.get('GB_LIST')
                elif 'Accession' in row:
                    gene_bank_accession = row.get('Accession')
                else:
                    print('Key not found!')

                platform_data[gene_id] = gene_bank_accession.strip()

    print('Platform loaded:', len(platform_data))
    return platform_data


def read_data(fr, fw, gse_accession, platform, platform_dir):
    print('gse_accession:', gse_accession, ', platform:', platform)
    platform_data = read_platform_data(platform, platform_dir)

    csv_reader = csv.reader(fr, delimiter='\t')
    csv_writer = csv.writer(fw, delimiter='\t')
    for row in csv_reader:
        if not row:
            continue # Skip emtpy lines

        if row[0] == 'ID_REF':
            row.insert(0, 'GENE')
            csv_writer.writerow(row) # write header
            continue

        # print(row)
        gene_id = row[0]
        gene_id_converted = platform_data[gene_id]
        # print('ID=', gene_id, 'New ID= ', gene_id_converted)
        tokens = gene_id_converted.split(' ')
        if len(tokens) > 1:
            print('tokens=', tokens)
        row.insert(0, gene_id_converted)
        csv_writer.writerow(row)

    print('File harmonized written')
    pass


def main(args):

    filename = os.path.join(args.dir, 'series_and_platforms.csv')
    series_and_platforms = read_series_and_platforms(filename)

    for entry in os.scandir(args.dir):
        if not entry.is_file():
            continue

        filename = entry.name
        # print(filename)

        if filename.endswith('.csv') or filename.endswith('_harmonized'):
            continue

        path = entry.path
        path_write = entry.path + '_harmonized'
        with open(path, 'r') as fr, open(path_write, 'w') as fw:
            tokens = filename.split('_')
            gse_accession = tokens[0]
            # print(gse_accession)

            read_data(fr, fw, gse_accession, series_and_platforms[gse_accession], args.platform_dir)

        # gse = geo.get_GEO_file('GSE37645')

        # full_filename = os.path.join(args.dir, filename)
        # gse = geo.parse_GSE(filename, {})
        # print(gse)
        # print(gse.get_accession())
        # print(gse.show_metadata())


if __name__ == '__main__':
    # Configure logging
    log_name = 'processing_file_names.log'
    logging.basicConfig(filename=log_name, filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s   %(message)s')
    logging.getLogger().addHandler(logging.StreamHandler())

    parser = argparse.ArgumentParser(description='Harmonize gene names among Series Matrix Files.')
    parser.add_argument('--dir', help='Directory', required=True)
    parser.add_argument('--platform_dir', help='Directory', required=True)

    args = parser.parse_args()
    main(args)