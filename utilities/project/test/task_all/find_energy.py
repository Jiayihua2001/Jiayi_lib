#!/usr/bin/env python
# Prompt the user for the file name
file_name = input("Enter the name of the file: ")

# Prompt the user for the target sentence
target_sentence = "Final output of selected total energy values:"

# Prompt the user for the number of lines to extract before and after the target sentence
context_lines =int(input("enter the lines you want to read,(usually 12)"))

try:
    # Open the file in read mode
    with open(file_name, 'r') as file:
        lines = file.readlines()

        found = False
        for line_number, line in enumerate(lines, start=-1):
            if target_sentence in line:
                found = True
                start_index = max(0, line_number-1)
                end_index = min(len(lines), line_number + context_lines + 1)

                # Extract and print the surrounding context
                context = lines[start_index:end_index]
                print(f"Found target sentence at line {line_number}:\n")
                print("".join(context))

                break

        if not found:
            print("Target sentence not found in the file.")

except FileNotFoundError:
    print(f"File '{file_name}' not found.")
except Exception as e:
    print(f"An error occurred: {str(e)}")

