#!/usr/bin/env python3

import sys
import tomllib

if len(sys.argv) == 3:
    pars = sys.argv[2]
    with open(sys.argv[1], "rb") as f:
        data = tomllib.load(f)
elif len(sys.argv) == 2:
    pars = sys.argv[1]
    with open("holi_params.toml", "rb") as f:
        data = tomllib.load(f)
else:
    print("Usage: get_params.py <path_to_toml_file> <key>")
    print("Usage: get_params.py  <key>")
    sys.exit(1)

# print(data)
# print(pars)
# print(pars.split("."))
value = data
for keys in pars.split("."):
    value = value[keys]

print(value)