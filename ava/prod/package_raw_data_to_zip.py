"""
Package HDX raw data into zip files, threaded, then checking if the zip files contain all files from the original folders.

thx Claude & ChatGPT.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue
from threading import Event
import zipfile
from tqdm.auto import tqdm
from pathlib import Path

# %%

src_dir = Path(r"C:\Users\jhsmi\_TEMP_SPLIT_DATA")
folders = [d for d in src_dir.iterdir() if d.is_dir()]

folder_tasks = {
    raw_f: [f for f in raw_f.glob("**/*") if f.is_file()] for raw_f in folders
}

# %%


def worker(task_id: int, folder: Path, q: Queue, stop: Event):
    """Do some work and report progress to the queue."""
    files = folder_tasks[folder]
    zip_file_path = folder.parent / (folder.name + ".zip")

    try:
        with zipfile.ZipFile(
            zip_file_path, "w", compression=zipfile.ZIP_DEFLATED
        ) as zipf:
            # zip a file
            for file in tqdm(files):
                zipf.write(file, file.relative_to(folder))
                if stop.is_set():
                    return  # early exit if something failed elsewhere
                q.put((task_id, 1))  # report +1 step
    except Exception as e:
        q.put((task_id, e))  # send the exception up
        raise


def run():
    num_tasks = len(folder_tasks)
    q = Queue()
    stop = Event()

    # Create one progress bar per task
    bars = []
    for i, (raw_f, files) in enumerate(folder_tasks.items()):
        bars.append(
            tqdm(
                total=len(files),
                position=i,  # stack bars vertically
                leave=True,  # keep bars after completion
                desc=f"Task {i + 1}",
            )
        )

    futures = []
    with ThreadPoolExecutor(max_workers=num_tasks) as ex:
        for i, (raw_f, files) in enumerate(folder_tasks.items()):
            futures.append(ex.submit(worker, i, raw_f, q, stop))

        finished = 0
        done_flags = [0] * num_tasks
        # Keep consuming progress until all bars reach total
        while finished < num_tasks:
            item = q.get()
            task_id, payload = item

            # If a worker sent an exception, stop everyone and re-raise later
            if isinstance(payload, Exception):
                stop.set()
                # Drain queue to avoid deadlocks
                while not q.empty():
                    q.get_nowait()
                break

            # Normal incremental update
            bars[task_id].update(payload)

            # Track completion when a bar reaches its total
            if (not done_flags[task_id]) and (bars[task_id].n >= bars[task_id].total):
                done_flags[task_id] = 1
                finished += 1

    # Surface any exceptions from threads
    for f in as_completed(futures):
        f.result()  # will re-raise if worker failed

    for b in bars:
        b.close()


# %%

run()

# %%

import os


def verify_zip_against_folder(zip_path: str, folder_path: str) -> tuple[bool, dict]:
    """
    Verify that a zip file contains all files from the original folder.

    Args:
        zip_path (str): Path to the zip file
        folder_path (str): Path to the original folder

    Returns:
        Tuple[bool, dict]: (is_complete, report_dict)
        - is_complete: True if zip contains all files from folder
        - report_dict: Dictionary with detailed comparison results
    """

    # Initialize report
    report = {
        "zip_exists": False,
        "folder_exists": False,
        "folder_files": [],
        "zip_files": [],
        "missing_from_zip": [],
        "extra_in_zip": [],
        "file_count_folder": 0,
        "file_count_zip": 0,
        "is_complete": False,
    }

    # Check if zip file exists
    if not os.path.exists(zip_path):
        report["error"] = f"Zip file not found: {zip_path}"
        return False, report
    report["zip_exists"] = True

    # Check if folder exists
    if not os.path.exists(folder_path):
        report["error"] = f"Folder not found: {folder_path}"
        return False, report
    report["folder_exists"] = True

    try:
        # Get files from folder (only files, not subdirectories, at root level)
        folder_files = set()
        for item in os.listdir(folder_path):
            item_path = os.path.join(folder_path, item)
            if os.path.isfile(item_path):
                folder_files.add(item)

        report["folder_files"] = sorted(list(folder_files))
        report["file_count_folder"] = len(folder_files)

        # Get files from zip (inside the folder with same name)
        folder_name = os.path.basename(folder_path)
        zip_files = set()
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            for zip_info in zip_ref.filelist:
                # Look for files inside the folder with the same name as the original folder
                if (
                    zip_info.filename.startswith(f"{folder_name}/")
                    and not zip_info.is_dir()
                ):
                    # Extract just the filename (remove the folder prefix)
                    file_name = zip_info.filename[len(folder_name) + 1 :]
                    # Only include files directly in that folder (no subdirectories)
                    if "/" not in file_name and file_name:
                        zip_files.add(file_name)

        report["zip_files"] = sorted(list(zip_files))
        report["file_count_zip"] = len(zip_files)

        # Compare file sets
        missing_from_zip = folder_files - zip_files
        extra_in_zip = zip_files - folder_files

        report["missing_from_zip"] = sorted(list(missing_from_zip))
        report["extra_in_zip"] = sorted(list(extra_in_zip))

        # Determine if complete
        is_complete = len(missing_from_zip) == 0
        report["is_complete"] = is_complete

        return is_complete, report

    except zipfile.BadZipFile:
        report["error"] = f"Invalid or corrupted zip file: {zip_path}"
        return False, report
    except Exception as e:
        report["error"] = f"Error during verification: {str(e)}"
        return False, report


def print_verification_report(report: dict, zip_path: str, folder_path: str):
    """
    Print a human-readable verification report.

    Args:
        report (dict): Report dictionary from verify_zip_against_folder
        zip_path (str): Path to zip file
        folder_path (str): Path to folder
    """
    print(f"\n=== ZIP VERIFICATION REPORT ===")
    print(f"Zip file: {zip_path}")
    print(f"Folder: {folder_path}")
    print(f"")

    if "error" in report:
        print(f"❌ ERROR: {report['error']}")
        return

    print(f"Files in folder: {report['file_count_folder']}")
    print(f"Files in zip: {report['file_count_zip']}")
    print(f"")

    if report["is_complete"]:
        print("✅ ZIP IS COMPLETE - All files are included!")
    else:
        print("❌ ZIP IS INCOMPLETE")

        if report["missing_from_zip"]:
            print(f"\n📁 Missing from zip ({len(report['missing_from_zip'])} files):")
            for file in report["missing_from_zip"]:
                print(f"  - {file}")

        if report["extra_in_zip"]:
            print(f"\n📦 Extra files in zip ({len(report['extra_in_zip'])} files):")
            for file in report["extra_in_zip"]:
                print(f"  + {file}")


# %%

folders
zips = [f.parent / (f.name + ".zip") for f in folders]

# %%

f, z = folders[0], zips[0]

for f, z in zip(folders, zips):
    print(f"\nVerifying {z} against {f}")
    is_complete, report = verify_zip_against_folder(str(z), str(f))
    print_verification_report(report, str(z), str(f))


# %%
