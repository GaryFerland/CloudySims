import sys
import os
import subprocess
import glob
import requests
import tarfile
import platform
import asyncio
import shutil
import subprocess
import re
from datetime import datetime

"""
Gold To Do Version: Making a Cloudy Release step by step instructions

Required packages: doxygen, pyppeteer, pdflatex (script will run through this as well)

    1. Do a merge from master to release branch. Two options to do this:

    Option 1: Do it manually
        1.1 get a fresh copy of master
                >> git pull
                >> git checkout master
                >> git fetch origin

        1.2 Update the copy right year
                >> find ./ -type f -exec sed -i -e 's/1978-2023/1978-2025/g' {} \;
            validate the changes
                >> grep -rnw . -e '1978' | grep -v 2025 | grep -v Percival | grep -v Draine

        1.3 Merge needed branches to master and make sure these changes are described in
           NewCXX and the release paper.

        1.4 Switch to the release branch
                >> git switch release

        1.5 do a fast forward merge from master onto the release branch with preference to
           incoming changes if there are any merge conflicts.

                >> git merge master --ff-only
            If merge fails with conflicts:
                >> git merge master
                >> git merge -X theirs master
                >> git push origin release
        The merge_master_to_release.sh does the steps till here.

    Option 2: Start on master and run the following shell script which goes through the above steps:
            >> chmod +x merge_master_to_release.sh
            >> ./merge_master_to_release.sh

            Note: merge_master_to_release.sh:31 does the push, this needs to be tested or commented out.

    Bring up two terminals. You will need one terminal to run this release preparer script,
    and another to check outputs, and follow the instructions provided.

    2. Run the release script from the root of the release branch

            >> python scripts/cloudy_release_preparer.py
        The script will ask if tsuite has been run, enter "n" to run the tsuite. Then
        once the tsuite is has run clean, come back, re-run the script and enter "y".

    3. Clean tsuite and source:
            >> cd source
            >> make dist clean

            >> cd ../tsuite
            >> ./clean_tsuite.pl

    4. Commit the changes

    5. Copy Doxygen to the data area
       NOTE: These should not go in the tarball

       First make a new subdirectory under doxygen for the new release
            >> ssh <user-name>@nublado.org
            >> cd /var/www/webapps/data_area/doxygen/
            >> mkdir c25.00
       Exit nublado.org
            >> [CTRL d]
       Copy the doxygen tree to the new subdirectory created in nublado
            >> rsync -a doxygen/html/ <user-name>@nublado.org:/var/www/webapps/data_area/doxygen/c25.00/

    6. Copy the release tarball to nublado
       (this script creates one automatically once all directories have been prepped sucessfully)
            >> rsync -avz c25.00.tar.gz <user-name>@nublado.org:/var/www/webapps/data_area/cloudy_releases/c25/

    7. Tag the latest release branch commit
"""

def check_packages():
    # Make sure Doxygen is installed
    # This checks if doxygen is installed, by running a version check
    try:
        result = subprocess.run(['doxygen', '--version'])
    except:
        print(f"Doxygen is not installed or not found in the system path.")

        install_success = False
        print("Attempting to download doxygen...")
        if platform.system() == "Darwin": # macOS
            try:
                subprocess.check_call(['brew', 'install', 'doxygen'])
                print("Doxygen installed!")
                install_success = True
            except:
                print("Failed to install Doxygen using Homebrew.")
                install_success = False

        # If error returned on doxygen version check, attempt to download doxygen
        if platform.system() != "Darwin" or install_success == False:
            try:
                subprocess.check_call(['sudo', 'apt-get', 'update'])
                subprocess.check_call(['sudo', 'apt-get', 'install', '-y', 'doxygen'])
                print("Doxygen installed!")
                install_success = True
            except:
                print("Failed to install Doxygen using sudo apt-get.")
                print("Please install doxygen manually before continuing. Aborting script.")
                return -1

    # pdfkit needed later for converting .htm file to .pdf file
    try:
        from pyppeteer import launch
    except:
        print(f"pyppeteer is not installed.")

        install_success = False
        print("Attempting to install pyppeteer...")
        try:
            subprocess.check_call(['pip', 'install', 'pyppeteer'])
            print("pyppeteer installed!")
        except:
            print("Failed to install pyppeteer using pip.")
            print("Please install pyppeteer manually before continuing. Aborting script.")
            return -1

    # Check if pdflatex is installed by doing a version check
    try:
        subprocess.run(["pdflatex", "--version"])
    except:
        print(f"pdflatex is not installed.")

        print("Attempting to install pdflatex...")
        if platform.system() == "Darwin": # macOS
            try:
                subprocess.run(["brew", "install", "--cask", "mactex"])
                print("mactex installed! Please re-start terminal and re-start script.")
            except:
                print("Failed to install mactex using Homebrew.")
                print("Please install mactex manually before continuing. Aborting script.")
                return -1

        elif platform.system() != "Darwin":
            try:
                subprocess.run(["sudo", "apt", "install", "texlive"])
                print("texlive installed! Please re-start terminal and re-start script.")
            except:
                print("Failed to install texlive using sudo apt-get.")
                print("Please install texlive manually before continuing. Aborting script.")
                return -1

    return 0

def update_copyright_year():
    # Check Copyright year and update to current year
    pattern = re.compile(r"1978-20\d{2}")
    current_year = str(datetime.now().year)
    replacement = f"1978-{current_year}"

    folder_path = "source/"
    for root, dirs, files in os.walk(folder_path):
        for filename in files:
            filepath = os.path.join(root, filename)
            try:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()

                # Change copyright year only if it needs to
                matches = pattern.findall(content)
                if matches and any(match != replacement for match in matches):
                    new_content = pattern.sub(replacement, content)
                    with open(filepath, 'w', encoding='utf-8') as f:
                        f.write(new_content)
            except Exception as e:
                print(f"Skipping {filepath}: {e}")

    # Also make this replacement in docs/latex/frontis_common.tex
    pattern = re.compile(r"1978--20\d{2}")
    replacement = f"1978--{current_year}"
    filepath = "docs/latex/frontis_common.tex"
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    matches = pattern.findall(content)
    if matches and any(match != replacement for match in matches):
        new_content = pattern.sub(replacement, content)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)

def update_readme(cloudy_release):
    if "_" in cloudy_release:
        print("Please review and update README.md in the root directory.")
        with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
            f.write("readme\n")
        return

    readme_file = "README.md"
    new_version = cloudy_release[1:]
    if "00" in cloudy_release:
        new_version = new_version.split(".")[0]
    year_suffix = str(datetime.now().year)[-2:]

    replacements = {
        r"The current version of Cloudy is C\d+, released in \d{4}.": f"The current version of Cloudy is C{new_version}, released in {datetime.now().year}.",
        r"\[here\]\(https://gitlab\.nublado\.org/cloudy/cloudy/-/wikis/NewC\d+\)": f"[here](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/NewC{year_suffix})"
    }

    with open(readme_file, 'r', encoding="utf-8") as f:
        lines = f.read()
    for pattern, replacement in replacements.items():
        lines = re.sub(pattern, replacement, lines)
    with open(readme_file, "w", encoding="utf-8") as f:
        f.write(lines)

    print("Please update the bib link to the latest release paper in README.md")
    readme_update_success = input(" Is the README.md file in the root directory up-to-date (y/n)? ")

    if readme_update_success == "y":
        with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
            f.write("readme\n")
    else:
        print("README.md not ready for release.")
        return


def prep_source(cloudy_release):
    os.chdir("./source/")

    current_dir = os.getcwd()
    print("Entered", current_dir)

    # This checks on any header files that are not used.
    command_args = ["./list_headers.pl"]
    print(f"\n Running source{command_args[0][1:]}")
    subprocess.run(command_args)

    header_summary = glob.glob(f"{current_dir}/headers.txt")[0]
    header_file_list = glob.glob(f"{current_dir}/listfiles.list")[0]
    print(f" Written to \n{header_summary} \n{header_file_list}")

    # This prunes duplicate headers from source files.
    command_args = ["./uninclude-headers.pl"]
    print(f"\n Pruning delicate headers from source files with {command_args[0]}")
    subprocess.run(command_args)

    # This extracts all atomic data references.
    command_args = ["./doc_atomic_data.pl"]
    print("\n Running ", command_args[0])
    subprocess.run(command_args)

    # Update the CLD_MAJOR, CLD_MINOR, CLD_BETA in version.cpp
    rc = None if len(cloudy_release.split("_")) == 1 else cloudy_release.split("_")[-1]
    CXX = cloudy_release.split(".")[0]
    CXX = CXX.lower()
    rev = cloudy_release.split(".")[1]
    if rc: rev = rev.split("_")[0]
    new_major = int(CXX.split("c")[-1])
    new_minor = int(rev)
    new_beta = 1 if rc else 0
    replacements = {
        r"(static\s+const\s+int\s+CLD_MAJOR\s*=\s*)\d+;": f"static const int CLD_MAJOR = {new_major};",
        r"(static\s+const\s+int\s+CLD_MINOR\s*=\s*)\d+;": f"static const int CLD_MINOR = {new_minor};",
        r"(static\s+const\s+int\s+CLD_BETA\s*=\s*)\d+;": f"static const int CLD_BETA = {new_beta};"
    }
    version_file = "version.cpp"
    # Read in version.cpp
    with open(version_file, "r", encoding="utf-8") as f:
        content = f.read()
    # Make the replacement with the new values for CLD_MAJOR, CLD_MINOR, CLD_BETA
    for pattern, replacement in replacements.items():
        content = re.sub(pattern, replacement, content)
    with open(version_file, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Version numbers in {version_file} updated to")
    print(f" CLD_MAJOR={new_major}, CLD_MINOR={new_minor}, CLD_BETA={new_beta}.")

    # Update the release date to current date in date.h
    date_file = "date.h"
    base_year = 2000
    today = datetime.now()
    year = today.year - base_year + 100 # So 2025 -> 125
    month = today.month - 1        # January is 0
    day = today.day
    replacement_dates = {
        r'#define\s+YEAR\t\d+':  f'#define YEAR\t{year}',
        r'#define\s+MONTH\t\d+': f'#define\tMONTH\t{month}',
        r'#define\s+DAY\t\d+':   f'#define\tDAY\t{day}'
    }
    with open(date_file, 'r', encoding="utf-8") as f:
        lines = f.read()
    for pattern, replacement in replacement_dates.items():
        lines = re.sub(pattern, replacement, lines)
    with open(date_file, "w", encoding="utf-8") as f:
        f.write(lines)
    print(f"Date in {date_file} updated to")
    print(f" YEAR={year}, MONTH={month}, DATE={day}.")

    print("\nSource directory ready for release.\n")
    os.chdir("../")
    with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
        f.write("source\n")


def prep_doxygen(cloudy_release):
    os.chdir("./doxygen/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    doxygen_html = glob.glob(f"{current_dir}/html/index.html")
    if doxygen_html:
        print("\n Doxygen files found: ", "/".join(doxygen_html[0].split("/")[-3:]))
        run_doxygen = input("\n Run Doxygen (y/n)? ")
    else:
        run_doxygen = "y"

    if run_doxygen.lower() == "y":
        # This creates the Doxygen documentation
        command_args = ["doxygen", "Doxyfile"]
        print("\n Running ", command_args[0], command_args[1])
        subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        doxygen_html = glob.glob(f"{current_dir}/html/index.html")

    # If documentation .html file exists open for viewing
    if f"{current_dir}/html/index.html" in doxygen_html:
        print("Doxygen successfully configured.\n Openning doxygen html\n")
        subprocess.run(["open", doxygen_html[0]])
    else:
        print("Oops! Could not configure doxygen. Please run the following command manually.")
        print(" >> doxygen Doxyfile")
        return

    print("Next you need to copy the doxygen tree to Cloudy's data server.")
    cloudy_data_server = f"[nublado_username]@nublado.org:/var/www/webapps/data_area/doxygen/{cloudy_release}/"
    print(f"rsync -avz html/ {cloudy_data_server}")
    doxygen_complete = input("Has the doxygen tree been copied to nublado (y/n)? ")

    if doxygen_complete == "y":
        print("\nDoxygen directory ready for release.\n")
        os.chdir("../")
        with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
            f.write("doxygen\n")


def prep_data():
    os.chdir("./data/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    readme_data_file = "README_data.md"
    print(f"\n Please review and update the data/{readme_data_file}.")
    readme_edit_success = input(" Is docs/README_data.md up-to-date (y/n)? ")
    if readme_edit_success == "n":
        print("Aborting release prep script. Come back once data/README_data.md has been updated.")
        return

    # This asks user to make sure all compiled data files are up to date.
    command_args = ["./make_data.sh"]
    print(f"\n Running {command_args[0]}, to make sure compiled data files are up to date.")
    subprocess.run(command_args)

    # This asks user to make sure Cloudy citations are up-to-date
    print("\n Please review and update data/citation_cloudy.txt")
    citation_update_success = input(" Is data/citation_cloudy.txt up-to-date (y/n)? ")
    if citation_update_success == "n":
        print("Aborting release prep script. Comeb back once data/citation_cloudy.txt has been updated.")
        return
    else:
        # Then test it by creating and running a test with print citation command
        with open("citation_test.in", 'w') as file:
            file.write("test\n")
            file.write("print citation")
        command_args = ["../source/cloudy.exe", "-r", "citation_test"]
        try:
            subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("\n I created and ran a test for the \'print citation\' command.")
            # Ask user to check the output from the test
            print(" Please review data/citation_test.out")
            citation_test_success = input(" Enter \'continue\' if test looks good, otherwise enter \'error\' to abort: ")
            if citation_test_success == "error":
                print("Error: aborting script.")
                return
            os.remove("citation_test.in")
            os.remove("citation_test.out")
        except:
            print("Error: aborting, something went wrong running Cloudy executable.")
            return

    # Ask user to update the citations for the databases used by Cloudy.
    print("\n Please review and update data/citation_data.txt, the file needs to be updated with the latest database versions.")
    citationdata_update_success = input(" Is data/citation_data.txt up-to-date (y/n)? ")
    if citationdata_update_success == "n":
        print("Aborting release prep script. Comeb back once data/citation_data.txt has been updated.")
        return -1

    # This makes sure checksums.dat is up to date.
    # If you build Cloudy in one of the sys_xxxx directories you must temporarily
    # copy (or symlink) vh128sum.exe into source.
    vh128sum_executable = glob.glob("../source/vh128sum.exe")
    if vh128sum_executable:
        command_args = ["../scripts/generate_checksums.sh"]
        print(f"\n Running {command_args[0]} to update checksums.dat")
        subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        print("A smoke test will now start to make sure checksums ran ok.")
        run_smoketest = input("Run smoke test (y/n) ?")
        if run_smoketest == "y":
            command_args = ["../source/cloudy.exe", "-e", "test"]
            subprocess.run(command_args)
    else:
        print("Could not find /source/vh128sum.exe. If you build Cloudy in one of the sys_xxxx")
        print("directories, you must temporarily copy (or symlink) vh128sum.exe into source.")
        print("Moving onto next directory.")
        return

    smoketest_success = input("Smoke test looks good (y/n) ?")

    if smoketest_success == "y":
        print("\nData directory ready for release.\n")
        os.chdir("../")
        with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
            f.write("data\n")
    else:
        print("Something went wrong. Data directory not ready for release.")
        return


async def convert_html_to_pdf(in_htm, out_pdf):
    browser = await launch(headless=True)
    page = await browser.newPage()
    await page.goto(in_htm, {'waitUntil': 'networkidle2'})  # local HTML file
    await page.pdf({'path': out_pdf})
    await page.pdf({
        'path': out_pdf,
        'margin': {
            'top': '1in',
            'right': '1in',
            'bottom': '1in',
            'left': '1in'
        }
        })
    await browser.close()


def prep_tsuite():
    os.chdir("./tsuite/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # Prepping tsuite/auto
    os.chdir("./auto/")
    print("\n Entered tsuite/auto/")

    # This creates a list of all test cases, including the input commands and a description of its purpose.
    command_args = ["./doc_tsuite.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, to update list of tests in tsuite/auto.")
    subprocess.run(command_args)

    # doc_tsuite.htm is output from above: contains a tab-delimited list of the files.
    # Open it for user to review
    subprocess.run(["open", "./doc_tsuite.htm"])
    print(" Review opened doc_tsuite.htm for tsuite/auto.")
    doc_tsuite_success = input(" Enter \'continue\' if doc_suites look good, otherwise enter \'error\' to abort: ")
    if doc_tsuite_success == "error":
        print("Error encountered, aborting release prep script.")
        return

    # This diagnoses results of tsuite run
    command_args = ["./checkall.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, to diagnose results of tsuite run.")
    files_to_check = ["wvlng.txt", "warnings.txt", "close.txt", "debug.txt", "serious.txt", "crashed.txt"]
    not_empty_files = []
    for file in files_to_check:
        if os.path.exists(file):
            if os.path.getsize(file) != 0:
                not_empty_files.append(file)
        else:
            print(f"Warning! {file} not found.")
    if not_empty_files:
        print("Following files are not empty. Please resolve and come back. Moving onto next directory.")
        return

    #The first character should be a sharp sign. This script lists all files that do not start with "#".
    # This is an error, and may indicate that the header was not properly produced.
    command_args = ["./CheckPunchSharp.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]} to make sure save files start with a header saying what the column indicates.")
    subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Prepping tsuite/slow with same steps as above
    os.chdir("../slow/")
    print("\n Entered tsuite/slow/")
    command_args = ["./doc_tsuite.pl"]
    print(f"\n Running tsuite/slow/{command_args[0][2:]}, to update list of tests in tsuite/slow.")
    subprocess.run(command_args)

    subprocess.run(["open", "./doc_tsuite.htm"])
    print("Review opened doc_tsuite.htm for tsuite/slow.")
    doc_tsuite_success = input(" Enter \'continue\' if doc_suites look good, otherwise enter \'error\' to abort: ")
    if doc_tsuite_success == "error":
        print("Error encountered, aborting release prep script.")
        return

    command_args = ["./checkall.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, to diagnose results.")
    files_to_check = ["wvlng.txt", "warnings.txt", "close.txt", "debug.txt", "serious.txt", "crashed.txt"]
    not_empty_files = []
    for file in files_to_check:
        if os.path.exists(file):
            if os.path.getsize(file) != 0:
                not_empty_files.append(file)
        else:
            print(f"Warning! {file} not found.")
    if not_empty_files:
        print("Following files are not empty. Please resolve and come back. Moving onto next directory.")
        return

    # Note: tsuite/slow/ does not have CheckPunchSharp.pl

    os.chdir("../")

    doctsuite_files = ["auto/doc_tsuite.pdf", "slow/doc_tsuite.pdf"]

    doctsuite_pdfs = []
    for file in doctsuite_files:
        found_file = glob.glob(file)
        if found_file: doctsuite_pdfs.append(found_file)

    print("\n pdf files found in tsuite: ", doctsuite_pdfs)

    if len(doctsuite_pdfs) == 2:
        html_to_pdf = input("\n Convert doc_tsuite.htm files to pdfs to be included in Hazy2 (y/n)? ")
    else:
        html_to_pdf = "y"

    if html_to_pdf == "y":
        print("\n Creating pdfs from doc_tsuite.htm files...")
        asyncio.get_event_loop().run_until_complete(convert_html_to_pdf(f"file://{current_dir}/auto/doc_tsuite.htm", f"{current_dir}/auto/doc_tsuite.pdf"))
        asyncio.get_event_loop().run_until_complete(convert_html_to_pdf(f"file://{current_dir}/slow/doc_tsuite.htm", f"{current_dir}/slow/doc_tsuite.pdf"))

        auto_doc = glob.glob(doctsuite_files[0])
        slow_doc = glob.glob(doctsuite_files[1])
        if not auto_doc or not slow_doc:
            print(" Something went wrong. Could not find auto/doc_tsuite.pdf and slow/doc_tsuite.pdf.")
            return

    # TODO: Find coverage run, what script does this? 

    os.chdir("./programs/")
    print("\n Entered tsuite/programs/")

    # This will run the set of programs in run_programs.dat with either gcc or icc
    command_args = ["./run_programs.pl"]
    print(f"\n Running tsuite/programs/{command_args[0][2:]}")
    subprocess.run(command_args)
    # TODO: what is needed to be checked here? Add test.

    print("Tsuite directory ready for release.\n")
    os.chdir("../../")
    with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
        f.write("tsuite\n")


def prep_scripts():
    os.chdir("./scripts/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # This creates a list of all test cases, including the input commands and a description of its purpose.
    command_args = ["./prep-bib.sh"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, gathers atomic data references from the Stout database into LaTeX/PDF format.")
    subprocess.run(command_args)

    print("Scripts directory ready for release.\n")
    os.chdir("../")
    with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
        f.write("scripts\n")


def prep_docs():
    os.chdir("./docs/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # This makes sure LineLabels.txt and SpeciesLabels.txt are up-to-date
    linelable_input_script = "LineLabels"
    print(f"\n Running docs/{linelable_input_script}.in")
    subprocess.run(["../source/cloudy.exe", "-r", linelable_input_script])
    outfile = glob.glob(f"{linelable_input_script}.out")
    if not outfile:
        print(f"Aborting script. Something went wrong running {linelable_input_script}.in")
        return

    os.chdir("./latex/")
    print("\n Entered", "docs/latex/")

    # Ask user to update citation in Hazy
    print("\n Please review and update citation in item \'CloudyReview\' in docs/latex/common/bibliography2.bib.")
    citation_update_success = input(" Enter \'continue\' once citation is updated, otherwise enter \'error\' to abort: ")
    if citation_update_success == "error":
        print("Encoutered error, aborting script.")
        return -1
    
    # Update table atomic data sources.
    print("\n Please review and update table atomic data sources.")
    datasource_update_success = input(" Enter \'continue\' once table atomic data sources are updated, otherwise enter \'error\' to abort: ")
    if datasource_update_success == "error":
        print("Encoutered error, aborting script.")
        return -1

    print("\n Check whether in-press papers in references have appeared. and update Hazy.")
    datasource_update_success = input(" Enter \'continue\' once in-press papers are updated, otherwise enter \'error\' to abort: ")
    if datasource_update_success == "error":
        print("Encoutered error, aborting script.")
        return -1
    
    print("\n Checking for remaining TODOs' in latex")
    result = subprocess.run(["grep", "-r", "--include=\"*.tex\"", "TODO", "."], capture_output=True, text=True)
    # Print the result
    if result.returncode == 0:
        print("Matches found: please update TODOs")
        todos_update = input(" Enter \'continue\' once TODOs are updated, otherwise enter \'error\' to abort: ")
        print(result.stdout)
    elif result.returncode == 1:
        print("No TODOs found. Moving on...")
    else:
        print("Encoutered error, aborting script.")
        return -1

    pdf_files = ["./hazy1/hazy1.pdf", "./hazy2/hazy2.pdf", "./hazy3/hazy3.pdf", "./QuickStart/QuickStart.pdf"]

    hazy_pdfs = []
    for file in pdf_files:
       found_file = glob.glob(file)
       if found_file: hazy_pdfs.append(found_file)

    print("\n Hazy pdfs found in docs/latex/: ", hazy_pdfs)
    # This provides an option to skip re-compiliing hazy pdfs if they have already been
    #  compiled in a release script run previously.
    if len(hazy_pdfs) == 4:
        compile_hazy = input("\n Recompile pdfs (y/n)? ")
    else:
        compile_hazy = "y"

    if compile_hazy.lower() == "y":
        print("\n Compiling Hazy latex files.")
        # This creates pdfs in each hazy subfolder
        command_args = ["./CompileAll.pl"]
        print(f"\n Running docs/latex/{command_args[0][2:]}, to creating Hazy pdf files.")
        subprocess.run(command_args)

    # Copy hazy1.pdf, hazy2.pdf, hazy3.pdf, and QuickStart.pdf to top of docs directory
    for file in pdf_files:
        shutil.copy(file, f"../{file.split("/")[-1]}")

    print("Docs directory ready for release.\n")
    os.chdir("../../")
    with open("cloudy_file_prep_log.txt", 'a', encoding='utf-8') as f:
        f.write("docs\n")


def main():
    print("\nFirst make sure that the executable has recently been made...")
    os.chdir("./source/")
    num_cpus = os.cpu_count()
    command_args = ["make", "-j", f"{num_cpus}"]
    subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    os.chdir("../")

    print("\nBefore we get started, the full tsuite must be run.")
    tsuite_run = input("Full tsuite has been run? (y/n) Warning: entering \'n\' will start tsuite run. ")
    if tsuite_run.lower() == "n":
        os.chdir("./tsuite/")
        subprocess.run(["./run_parallel.pl"])
        os.chdir("../")
        return
    elif tsuite_run.lower() == "y":
        update_copyright_year()

        # Write release prep log
        log_file = glob.glob("./cloudy_file_prep_log.txt")
        if not log_file:
            with open("./cloudy_file_prep_log.txt", 'w', encoding='utf-8') as f:
                f.write("# Following directories are prepped for release:\n")
            print("Release prep log file created: cloudy_file_prep_log.txt")
        cloudy_release = input("Enter cloudy release version number (e.g. \'c25.00\'): ")
        with open("./cloudy_file_prep_log.txt", 'r', encoding='utf-8') as f:
            release_log = f.read()
        if "readme" not in release_log: update_readme(cloudy_release)
        if "source" not in release_log: prep_source(cloudy_release)
        if "doxygen" not in release_log: prep_doxygen(cloudy_release)
        if "data" not in release_log: prep_data()
        if "tsuite" not in release_log: prep_tsuite()
        #if "scripts" not in release_log: prep_scripts()
        if "docs" not in release_log: prep_docs()

        with open("./cloudy_file_prep_log.txt", 'r', encoding='utf-8') as f:
            release_log = f.read()
        if all(s in release_log for s in ["readme", "source", "doxygen", "data", "tsuite", "docs"]):
            print("All directories prepped.")

            # Define your parameters
            output_file = f"{cloudy_release}.tar.gz"
            branch_name = "release"
            prefix_dir = f"{cloudy_release}/"

            # Run the git archive command
            print(f"Now creating {cloudy_release} tarball...")
            subprocess.run([
                "git", "archive",
                "--format=tar.gz",
                f"--output={output_file}",
                f"--prefix={prefix_dir}",
                branch_name
            ], check=True)

            print(f"Tarball created: {output_file}")

            if os.path.isfile(log_file[0]):
                os.remove(log_file[0])
        else:
            print("I did not create a release tarball since some directories failed to be prepped.")
    else:
        print("Aborting release prep script!")
        return

package_success = check_packages()
if package_success < 0: sys.exit(1)

from pyppeteer import launch

if __name__ == "__main__":
    main()

# TODO: add in a read and update docs/latex/version.tex to the latest version number
# TODO: some of these steps need to be in master, others only in release. This needs to be disentangled.
