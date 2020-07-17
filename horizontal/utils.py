import ujson as json
import os

def create_path(path):
    if not os.path.exists(path):
        os.makedirs(path)


def to_json(out_fname, data):
    with open(out_fname, "w") as wopen:
        json.dump(data, wopen, indent=4)


def from_json(fname):
    with open(fname, "r") as ropen:
        data = json.load(ropen)
    return data
