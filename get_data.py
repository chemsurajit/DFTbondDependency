import os
import requests
import logging
import zipfile

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

def download_files(download_path, base_url, download_option):
    """
    Download zip files from figshare.
    """
    logging.info("Downloading files...")

    if not os.path.exists(download_path):
        os.makedirs(download_path)

    if not os.path.exists(os.path.join(download_path, ".gitignore")):
        with open(os.path.join(download_path, ".gitignore"), "w") as f:
            f.write("*")

    if download_option == "all":
        files = {
            "csvs.zip": 19947446,
            "outputs.zip": 19942154,
            "xyzfiles.zip": 19780570,
            "TZP.zip": 19786672
        }
    if download_option == "partial":
        files = {
            "xyzfiles.zip": 19780570,
            "TZP.zip": 19786672
        }

    for i, (file_name, file_id) in enumerate(files.items()):
        file_url = base_url + str(file_id)
        file_dest = os.path.join(dest, file_name)
        if not os.path.exists(file_dest):
            logging.info("Downloading files: %d/%d --> %s" % (i+1, len(files), file_url))
            download_file(file_url, file_dest)
        else:
            logging.info("File already exists: %d/%d --> %s" %(i+1, len(files), file_dest))

    return

def extract_zip(src, dest):
    """
    Extract one zip file.
    """
    with zipfile.ZipFile(src, "r") as zip:
        zip.extractall(path=dest)
    return

def extract_zips(download_path, extract_path, download_option):
    """
    Extract all the zip files.
    """

    if download_option == "all":
        files = ["csvs.zip",
                 "outputs.zip",
                 "xyzfiles.zip",
                 "TZP.zip"
        ]
    if download_option == "partial":
        files = [
            "xyzfiles.zip",
            "TZP.zip"
        ]

    logging.info("Extracting files...")
    src = download_path
    assert os.path.exists(src)

    if os.path.exists(extract_path):
        logging.info("Directory already exists: %s" % extract_path)
        return

    for zipf in files:
        output_file_path = os.path.join(extract_path, zipf.split(".")[0])
        zip_file_path = os.path.join(download_path, zipf)
        if os.path.exists(zip_file_path):
            logging.info("Extracting file to: %s" % output_file_path)
            extract_zip(zip_file_path, output_file_path)
        else:
            logging.warning("File %s not found" % zip_file_path)

def get_arguments():


if __name__ == "__main__":
    download_path = os.path.abspath("./downloads")
    base_url = sys.argv[1]
    extract_path = os.path.abspath("./")
    download_files(download_path, base_url)
    extract_zips(download_path, extract_path)
