import json
import os
import sys
import requests
import logging
import zipfile
import argparse


def download_file(url, dest):
    """
    Download files from the url given. In this case, it is DTU Data.
    """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return


def download_files(download_path, links_dict, zipfiles):
    """
    Download zip files from figshare.
    """
    logging.info("Downloading files...")

    if not os.path.exists(download_path):
        os.makedirs(download_path)

    if not os.path.exists(os.path.join(download_path, ".gitignore")):
        with open(os.path.join(download_path, ".gitignore"), "w") as f:
            f.write("*")

    for i, file_name in enumerate(zipfiles):
        file_url = links_dict[file_name]
        logging.info("file url: %s" % file_url)
        file_dest = os.path.join(download_path, file_name+".zip")
        if not os.path.exists(file_dest):
            logging.info("Downloading files: %d/%d --> %s" % (i+1, len(zipfiles), file_url))
            download_file(file_url, file_dest)
        else:
            logging.info("File already exists: %d/%d --> %s" %(i+1, len(zipfiles), file_dest))
    return

def extract_zip(src, dest):
    """
    Extract one zip file.
    """
    with zipfile.ZipFile(src, "r") as zip:
        zip.extractall(path=dest)
    return

def extract_zips(download_path, extract_path, zipfiles):
    """
    Extract all the zip files.
    """

    to_logs_dir = ["TZP", "DZP", "SZ"]

    logging.info("Extracting files...")
    assert os.path.exists(download_path)

    for zipf in zipfiles:
        if zipf in to_logs_dir:
            output_file_path = os.path.join(extract_path, "logs")
        else:
            output_file_path = extract_path
        #output_file_path = os.path.join(extract_path, zipf.split(".")[0])
        zip_file_path = os.path.join(download_path, zipf+".zip")
        if os.path.exists(zip_file_path):
            logging.info("Extracting file %s to: %s" % (zipf+".zip", output_file_path))
            extract_zip(zip_file_path, output_file_path)
        else:
            logging.warning("File %s not found" % zip_file_path)
    return


def get_arguments():
    """
    Function to parse arguments.
    """
    parser = argparse.ArgumentParser("File to download data from the database.")
    parser.add_argument(
        "-all", "--all",
        help="Option to download all data. If present, all data will be downloaded. \
             if not present, only the xyzfiles and logfiles will be downloaded.",
        action="store_true"
    )
    parser.add_argument(
        "-links", "--links",
        help="json file containing all the links for the data.",
        type=str,
        required=True
    )
    parser.add_argument(
        '-log', '--log',
        type=str,
        required=False,
        default="info",
        choices=["debug", "info", "warning", "error", "critical"],
        help="Provide logging level. Default is warning."
    )
    return parser.parse_args()


if __name__ == "__main__":
    # all directory will be relative to the CWD by default.
    # downloads to $CWD/zips
    # log files extracted to $CWD/logs/TZP $CWD/logs/DZP...
    # final reaction outputs containing ens, bond changes to $CWD/outputs/
    # other already calculated csvs to $CWD/csvs
    # xyzfiles to $CWD/xyzfiles
    download_path = os.path.abspath("./zips")
    args = get_arguments()
    # log settings
    log_level = args.log.upper()
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        level=log_level,
        datefmt="%H:%M:%S",
    )
    #
    link_js = args.links
    with open(link_js) as fp:
        links_dict = json.load(fp)
    if args.all:
        files = ["csvs", "outputs", "xyzfiles", "TZP"]
    else:
        files = ["xyzfiles", "TZP"]
    extract_path = os.path.abspath("./")
    download_files(download_path, links_dict, files)
    extract_zips(download_path, extract_path, files)
