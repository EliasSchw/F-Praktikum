def writeLatexMacros(macro_name, value, unit=None, filepath='Paper/macros.tex'):
    """
    Writes or overrides LaTeX macros in the specified file.
    """

    macro_content = f"\\newcommand{{\\{macro_name}}}{{{value:.2e}"
    if unit:
        macro_content += f"\\,{unit}"
    macro_content += "}}\n"

    # Read the existing file content
    try:
        with open(filepath, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        lines = []

    # Remove existing macro with the same name
    lines = [line for line in lines if not line.strip().startswith(f"\\newcommand{{\\{macro_name}}}")] 

    # Append the new macro
    lines.append(macro_content)

    # Write back to the file
    with open(filepath, 'w') as file:
        file.writelines(lines)





