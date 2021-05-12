#!/usr/bin/env python3

header_filename = 'include/cxxopts.hpp'
version_variable_prefix="CXXOPTS_"
version = {}

try:
    with open(header_filename) as header_file:
        for line in header_file:
            if version_variable_prefix + "_VERSION_MAJOR" in line:
                version["major"] = line.split()[-1]
            if version_variable_prefix + "_VERSION_MINOR" in line:
                version["minor"] = line.split()[-1]
            if version_variable_prefix + "_VERSION_PATCH" in line:
                version["patch"] = line.split()[-1]

            if len(version.keys()) == 3:
                break

    if len(version.keys()) != 3:
        print("ERROR: Unable to find version!")
        exit(-1)

    version_string = "{major}.{minor}.{patch}".format(**version)
    print("Current version:", version_string)
except FileNotFoundError:
    print("Header 'include/cxxopts.hpp' not found!")
    pass
    

new_version_string = input("Enter new version: ")

url = ("https://raw.githubusercontent.com/jarro2783/cxxopts/v{}/{}"
       .format(new_version_string, header_filename))

import urllib.request

print("Fetching updated header at:", url)
try:
    urllib.request.urlretrieve(url, header_filename)
except Exception as err:
    print("ERROR: Unable to download updated header!")
    print(err)
    exit(-2)

