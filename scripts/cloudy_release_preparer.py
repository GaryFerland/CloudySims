import os
import subprocess
import glob
import requests
import tarfile

import subprocess



def prep_source():
    os.chdir("./source/")

    current_dir = os.getcwd()
    print("Entered", current_dir)
    command_args = ["./list_headers.pl"]
    print("\n Running ", command_args[0])
    subprocess.run(command_args)
    header_summary = glob.glob(f"{current_dir}/headers.txt")[0]
    header_file_list = glob.glob(f"{current_dir}/listfiles.list")[0]
    print(f" Written to \n{header_summary} \n{header_file_list}")

    command_args = ["./uninclude-headers.pl"]
    print(f"\n Pruning delicate headers from source files with {command_args[0]}")
    subprocess.run(command_args)

    command_args = ["./doc_atomic_data.pl"]
    print("\n Running ", command_args[0])
    subprocess.run(command_args)

    cloudy_executable = glob.glob(f"{current_dir}/cloudy.exe")
    if f"{current_dir}/cloudy.exe" not in cloudy_executable:
        num_cpus = os.cpu_count()
        command_args = ["make", "-j", f"{num_cpus}"]
        print("Making Cloudy executable for later use.")
        subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    os.chdir("../")
    print("\nSource directory ready for release.\n")


def prep_doxygen(cloudy_release):
    os.chdir("./doxygen/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # Try running the 'doxygen --version' command to check if Doxygen is installed
    try:
        result = subprocess.run(['doxygen', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # If the returncode is 0, Doxygen is installed
        if result.returncode == 0:
            result = {"message": f"Doxygen is installed. Version: {result.stdout.strip()}", "x": 0}
        else:
            result = {"message": f"Doxygen is not installed or not found in the system path.\n Error: {result.stderr.strip()}", "x": -1}
    except:
        result = {"message": f"Doxygen is not installed or not found in the system path.", "x": -1}

    print(result["message"])
    
    # If error returned on doxygen version check, attempt to download doxygen
    if result["x"] == -1:
        print("Attempting to download doxygen...")
        try:
            # Install doxygen using Homebrew
            subprocess.check_call(['brew', 'install', 'doxygen'])
            print("Doxygen installed!")
            install_success = True
        except subprocess.CalledProcessError as e:
            print("Failed to install Doxygen using Homebrew.")
            install_success = False

        if install_success == False:
            try:
                subprocess.check_call(['sudo', 'apt-get', 'update'])
                subprocess.check_call(['sudo', 'apt-get', 'install', '-y', 'doxygen'])
                print("Doxygen installed!")
                install_success = True                
            except subprocess.CalledProcessError as e:
                print("Failed to install Doxygen using apt-get.")
                install_success = False
        
        if install_success == False:
            print("Please install doxygen manually before continuing.")
            return

    command_args = ["doxygen", "Doxyfile"]
    print("\n Running ", command_args[0], command_args[1])
    subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    doxygen_html = glob.glob(f"{current_dir}/html/index.html")
    doxygen_latx = glob.glob(f"{current_dir}/latex/")
    print(doxygen_html)
    if f"{current_dir}/html/index.html" in doxygen_html:
        print("Doxygen successfully configured.\n Openning doxygen html\n")
        subprocess.run(["open", doxygen_html[0]])

    # FOLLOWING BIT NEEDS TO BE REVISED: 
    # I think we need to be in nebula
    # make sure to first test the if the directory is there
#    print("Next we need to copy the doxygen tree to Cloudy's data server.")
#    nublado_username = input("Please enter (https://data.nublado.org/) username: ")
#    #nublado_password = input("password: ")
#    #print(nublado_username, nublado_password)
#    cloudy_data_server = f"{nublado_username}@nublado.org:/var/www/webapps/data_area/doxygen/{cloudy_release}/"
#    command_args = ["rsync", "-avz", "html/", cloudy_data_server]
#    try:
#        print(f"Copying doxygen directory try to @nublado.org:/var/www/webapps/data_area/doxygen/{cloudy_release}")
#        subprocess.run(command_args)
#    except:
#        print("Could not copy doxygen directory to @nublado.org:/var/www/webapps/data_area/doxygen/")

    os.chdir("../")
    print("\nDoxygen directory ready for release.\n")


def prep_data():
    os.chdir("./data/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    readme_data_file = "README_data.md"
    print(f"\n Please review and update the data/{readme_data_file}.")
    readme_edit_success = input(" Enter \'continue\' when finished review and update, or enter \'error\' to abort release prep: ")

    vh128sum_executable = glob.glob("../source/vh128sum.exe")

    if readme_edit_success.lower() == "error":
        print("Error encountered, aborting release prep script.")
        return
    elif readme_edit_success.lower() == "continue" and vh128sum_executable != []:
        command_args = ["../scripts/generate_checksums.sh"]
        #print(vh128sum_executable[0])
        print(f"\n Running {command_args[0]} to update checksums.dat")
        subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # TODO: add test for checksums run success 

    command_args = ["./make_data.sh"]
    print(f"\n Running {command_args[0]}, to make sure compiled data files are up to date.")
    subprocess.run(command_args)

    print("\n Please review and update data/citation_cloudy.txt")
    citation_update_success = input(" Enter \'continue\' once Cloudy citations have been updated, otherwise enter \'error\' to abort: ")
    if citation_update_success == "error":
        print("Error: aborting script.")
    else:
        with open("citation_test.in", 'w') as file:
            file.write("test\n")
            file.write("print citation")
        command_args = ["../source/cloudy.exe", "-r", "citation_test"]
        try:
            subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("\n I created and ran a test for the \'print citation\' command.")
            print(" Please review data/citation_test.out")
            citation_test_success = input(" Enter \'continue\' if test looks good, otherwise enter \'error\' to abort: ")
            if citation_test_success == "error":
                print("Error encountered, aborting release prep script.")
                return
            subprocess.run(["rm", "citation_test*"])
        except:
            print("Error: aborting, something went wrong running Cloudy executable.")
            return

    print("\n Please review and update data/citation_data.txt, the file needs to be updated with the latest database versions.")
    citationdata_update_success = input("Enter \'continue\' once Cloudy citations have been updated, otherwise enter \'error\' to abort: ")
    if citationdata_update_success == "error":
        print("Error encountered, aborting release prep script.")
        return
    else:
        os.chdir("../")
        print("\nData directory ready for release.\n")


def prep_tsuite():
    os.chdir("./tsuite/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    # TSUITE/AUTO
    os.chdir("./auto/")
    print("\n Entered tsuite/auto/")
    command_args = ["./doc_tsuite.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]}, to update list of tests in tsuite/auto.")
    subprocess.run(command_args)

    subprocess.run(["open", "./doc_tsuite.htm"])
    print(" Review opened doc_tsuite.htm for tsuite/auto.")
    doc_tsuite_success = input(" Enter \'continue\' if doc_suites look good, otherwise enter \'error\' to abort: ")
    if doc_tsuite_success == "error":
        print("Error encountered, aborting release prep script.")
        return

    command_args = ["./CheckPunchSharp.pl"]
    print(f"\n Running tsuite/auto/{command_args[0][2:]} to make sure save files start with a header saying what the column indicates.")
    #The first character should be a sharp sign. This script lists all files that do not start with "#". This is an error, and may indicate that the header was not properly produced.
    subprocess.run(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # TSUITE/SLOW
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

    # Note: CheckPunchSharp.pl in tsuite/slow
    os.chdir("../")

    # Following script has problems, TODO: fix
#    print("\n Creating pdf from new do doc_tsuite.htm files to include in Hazy2.")
#    # Install needed libraries; TODO: following needs to be handled better if subprocess.run does not work in any instance
#    try:
#        import pdfkit
#    except ImportError:
#        subprocess.run(["pip", "install", "pdfkit"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#        import pdfkit
#    
#    try:
#        subprocess.run(["wkhtmltopdf","--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#    except:
#            try:
#                subprocess.run(["brew", "install", "wkhtmltopdf"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#            except:
#                subprocess.run(["sudo", "apt-get", "install", "wkhtmltopdf"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#
#    pdfkit.from_file("./auto/doc_tsuite.htm", "./auto/doc_tsuite.pdf")
#    pdfkit.from_file("./slow/doc_tsuite.htm", "./slow/doc_tsuite.pdf")
#    print(" New doc_tsuite.pdf files created from doc_tsuite.htm files")
#    #TODO: these need to go to Hazy2

    # TODO: Find coverage run, what script does this? 

    os.chdir("./programs/")
    print("\n Entered tsuite/programs/")
    command_args = ["./run_programs.pl"]
    print(f"\n Running tsuite/programs/{command_args[0][2:]}")
#    subprocess.run(command_args)
    # TODO: what is needed to be checked here? Add test.
    # TODO: enter progress bar using time.sleep(0.1)

    os.chdir("../../")
    print("Tsuite directory ready for release.\n")


def prep_docs():
    os.chdir("./docs/")
    current_dir = os.getcwd()
    print("Entered", current_dir)

    linelable_input_script = "LineLables"
    print(f"\n Running docs/{linelable_input_script}.in")
    subprocess.run(["../source/cloudy.exe", "-r", linelable_input_script])

    os.chdir("./latex/")
    print("\n Entered", "docs/latex/")

    print("\n Please review and update citation in item \'CloudyReview\' in docs/latex/common/bibliography2.bib.")
    citation_update_success = input(" Enter \'continue\' once citation is updated, otherwise enter \'error\' to abort: ")
    if citation_update_success == "error":
        print("Encoutered error, aborting script.")
        return
    
    print("\n Please review and update table atomic data sources.")
    citation_update_success = input(" Enter \'continue\' once table atomic data sources are updated, otherwise enter \'error\' to abort: ")
    if citation_update_success == "error":
        print("Encoutered error, aborting script.")
        return
    
    # TODO: Check whether in-press papers in references have appeared.

    print("\n Checking for remaining TODOs' in latex")
    subprocess.run(["grep", "-r", "--include=\"*.tex\"", "TODO", "."])

    hazy_pdfs = glob.glob("hazy*.pdf")
    print(hazy_pdfs)
    if "hazy1.pdf" in hazy_pdfs and "hazy2.pdf" in hazy_pdfs and "hazy3.pdf" in hazy_pdfs:
        compile_hazy = input("\n Hazy pdfs found in docs/latex/. Recompile pdfs (y/n)? ")
    else:
        compile_hazy = "y"

    if compile_hazy.lower() == "y":
        print("\n Compiling Hazy latex files.")
        try:
            subprocess.run(["pdflatex", "--version"])
        except:
            print("Your system does not have pdflatex. Wait a moment, attempting to install pdflatex...")
            os_sys = input("What OS system are you using (e.g. mac, linux)? ")
            if "mac" in os_sys:
                try:
                    subprocess.run(["brew", "install", "--cask", "mactex"])
                except:
                    print("Unable to install pdflatex, please install and try again.\nAborting script!")
                    return
            if "linux" in os_sys:
                try:
                    subprocess.run(["sudo", "apt", "install", "texlive"])
                except:
                    print("Unable to install pdflatex, please install and try again.\nAborting script!")
                    return

        command_args = ["./CompileAll.pl"]
        print(f"\n Running docs/latex/{command_args[0][2:]}, to creating Hazy pdf files.")
        subprocess.run(command_args)

        # TODO: Replace any old pdf versions of hazy in the top directory with the current version and add them to the branch. 
        # copy hazy1.pdf, hazy2.pdf, hazy3.pdf, and QuickStart.pdf to top of docs directory and commit to release prep repo
        # The original names, as in "hazy1.pdf", must not be changed since cross references rely on them.

    os.chdir("../")
    print("Docs directory ready for release.\n")


def main():
    print("Before we get started, the full tsuite must be run.")
    tsuite_run = input("Full tsuite has been run? (y/n) Warning: entering \'n\' will start tsuite run. ")
    if tsuite_run.lower() == "n":
        os.chdir("./tsuite/")
        subprocess.run(["./run_parallel.pl"])
        os.chdir("../")
    elif tsuite_run.lower() == "y":
        cloudy_release = input("Enter cloudy release version number (e.g. \'c25.00\'):")
        prep_source()
        prep_doxygen(cloudy_release)
        prep_data()
        prep_tsuite()
        prep_docs()
    else:
        print("Aborting release prep script!")
        return

if __name__ == "__main__":
    main()

    #TODO:
    # i. copy doxygen tree to data area, do this in nebula probably
    # ii. make the returns 'return -1', otherwise if success 'return 0' 
    #     then give success message of each successful directory at very end, 
    #      or unsuccessful directories at the very end.
    # iii. once all directories are successfully prepped, the script should make the tarball
    # iv. commit the changes, then do a fresh checkout once all directories prepped: 
    #       - tarball then clean tsuite/.out
    #       - clean source directory
    #   NOTE: rather than doing a fresh checkout its best to have two copies of the release_prep branch, 
    #         one where you will run this script, and build all the needed files, and one clean on to be made 
    #         into a tarball. This way you wont accidentally lose files made by this script in the fresh check
    #         out. 
    #
    #   ONCE THE SCRIPT IS RUN DO THE FOLLOWING:
    # 
    # tarball: 
    # git archive --format=tar.gz --output=c25.00_rc1.tar.gz --prefix=c25.00_rc1/ c25_release_prep
    #
    # Add tarball to data_area:
    # rsync -avz c25.00.tar.gz <user-name>@nublado.org:/var/www/webapps/data_area/cloudy_releases/c25/
    #
    # Add doxygen to the data_area:
    # NOTE: These should not go in the tarball
    # ssh <user-name>@nublado.org 
    # cd /var/www/webapps/data_area/doxygen/
    # mkdir c25.00
    # [CTRL d]
    # rsync -a doxygen/html/ cmgu228@nublado.org:/var/www/webapps/data_area/doxygen/c25.00/