import os

def replace_text_in_file(file_path, old_text, new_text):
    try:
        # Read the file
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Replace occurrences of the old text
        updated_content = content.replace(old_text, new_text)

        # Write back the updated content
        with open(file_path, 'w', encoding='utf-8') as file:
            file.write(updated_content)
        print(f"Updated: {file_path}")

    except Exception as e:
        print(f"Failed to process {file_path}: {e}")


def replace_text_in_all_files(root_dir, old_text, new_text):
    for foldername, subfolders, filenames in os.walk(root_dir):
        for filename in filenames:
            file_path = os.path.join(foldername, filename)

            # Process only text-based files
            if filename.endswith(('.txt', '.md', '.html', '.csv', '.json', '.xml', '.py')):  # Add more extensions if needed
                replace_text_in_file(file_path, old_text, new_text)


if __name__ == "__main__":
    # Set root directory
    root_directory = input("Enter the directory to start searching: ")
    old_string = "1978-2023"
    new_string = "1998-2025"

    replace_text_in_all_files(root_directory, old_string, new_string)
    print("Replacement completed.")
