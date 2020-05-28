#!/bin/env python
# file name: index_gen.py
# creates indexes for the release links on the genome center website
# executed in checkMD5_fastq.sh

from pathlib import Path
import argparse

PARSER = argparse.ArgumentParser()
PARSER.add_argument('--directory', '-d', metavar='path/to/directory/to/index',
                    help='directory that needs index file, and sub dirs')
ARGS = PARSER.parse_args()

DIRECTORY = Path(ARGS.directory)

def sizeof_fmt(num):
    for unit in ['B','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return f'{round(num,2)}{unit}'
        num /= 1024.0
    return f'{num}Y'

def index_dir(directory):
    index = directory / 'index.html'
    try:
        index.unlink()
    except FileNotFoundError:
        pass

    files_dirs = [fd for fd in directory.glob('[!.]*')]
    dirs = [d_ for d_ in sorted(files_dirs) if d_.is_dir()]
    files = [f_ for f_ in sorted(files_dirs) if f_.is_file()]

    with index.open('w') as index_file:
        print("<html><body> <p><ul>", file=index_file)
        for f_ in files:
            print(f"<li><a href={f_.name}>{f_.name}</a> &nbsp; &nbsp; ({sizeof_fmt(f_.stat().st_size)})</li>", file=index_file)
            f_.chmod(0o664)
        for d_ in dirs:
            print(f"<li><a href={d_.name}>{d_.name}</a> &nbsp; &nbsp; ({sizeof_fmt(d_.stat().st_size)})</li>", file=index_file)
        print("</ul></p></body></html>", file=index_file)
    index.chmod(0o664)

    for d_ in dirs:
        print(f"moving on to {d_.name}")
        index_dir(d_)
        d_.chmod(0o775)

if __name__ == '__main__':
    index_dir(DIRECTORY)
