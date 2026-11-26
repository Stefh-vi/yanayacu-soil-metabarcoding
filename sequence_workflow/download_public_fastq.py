"""List or download the 36 public FASTQ files mapped to the current nine samples.

Default behavior is a manifest-only preflight. An explicit --download and a new
output directory are required for network transfer. ENA gzip files are verified
against their ENA size and MD5, not the canonical vendor gzip SHA-256.
"""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import json
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=ROOT/'data/raw_manifest/analysis_fastq_manifest.tsv')
    parser.add_argument('--marker', choices=['16S', 'ITS', 'both'], default='both')
    parser.add_argument('--sample', action='append', help='Optional exact sample ID; repeat to select multiple samples')
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--output', type=Path, help='New directory for explicitly requested downloads')
    args = parser.parse_args()
    with args.manifest.open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert len(rows) == 36 and len({(r['sample_id'], r['marker'], r['direction']) for r in rows}) == 36
    samples = {r['sample_id'] for r in rows}
    if args.sample and not set(args.sample) <= samples:
        parser.error('Unknown sample ID')
    rows = [r for r in rows if (args.marker == 'both' or r['marker'] == args.marker) and (not args.sample or r['sample_id'] in args.sample)]
    for row in rows:
        u = urllib.parse.urlparse(row['ena_fastq_url'])
        assert u.scheme == 'https' and u.hostname == 'ftp.sra.ebi.ac.uk'
        assert len(row['ena_fastq_md5']) == 32 and int(row['ena_fastq_bytes']) > 0
    summary = {'mode': 'manifest_preflight', 'files': len(rows), 'bytes': sum(int(r['ena_fastq_bytes']) for r in rows),
               'samples': sorted({r['sample_id'] for r in rows}), 'markers': sorted({r['marker'] for r in rows}),
               'raw_reads_downloaded': False, 'vendor_gzip_equivalence': 'not established by public archive checksums'}
    if not args.download:
        print(json.dumps(summary, indent=2)); return
    if not args.output:
        parser.error('--download requires --output pointing to a new directory')
    destination = args.output.resolve()
    if destination.exists():
        parser.error('Output already exists; choose a new directory. No files will be overwritten.')
    destination.mkdir(parents=True, exist_ok=False)
    receipts = []
    for row in rows:
        filename = Path(urllib.parse.urlparse(row['ena_fastq_url']).path).name
        path = destination/row['marker']/filename
        path.parent.mkdir(exist_ok=True)
        part = path.with_suffix(path.suffix+'.part')
        md5 = hashlib.md5(); sha = hashlib.sha256(); size = 0
        with urllib.request.urlopen(row['ena_fastq_url'], timeout=90) as response, part.open('xb') as stream:
            while True:
                block = response.read(1024*1024)
                if not block: break
                stream.write(block); md5.update(block); sha.update(block); size += len(block)
        if size != int(row['ena_fastq_bytes']) or md5.hexdigest() != row['ena_fastq_md5']:
            raise RuntimeError('Archive integrity mismatch; incomplete file retained at '+str(part))
        part.rename(path)
        receipts.append({**row, 'downloaded_path': str(path.relative_to(destination)), 'bytes_verified': size,
                         'md5_verified': md5.hexdigest(), 'download_sha256': sha.hexdigest()})
        (destination/'DOWNLOAD_RECEIPT.json').write_text(json.dumps({'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'status': 'COMPLETE' if len(receipts) == len(rows) else 'IN_PROGRESS', 'files': receipts}, indent=2)+'\n', encoding='utf-8')
    summary.update(mode='download_completed', raw_reads_downloaded=True, files_verified=len(receipts))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
