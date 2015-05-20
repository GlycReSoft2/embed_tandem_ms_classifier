import argparse
import os
import sys
import subprocess

try:
    from urllib.request import urlopen
    from urllib.error import URLError, HTTPError
except ImportError:
    from urllib2 import urlopen, URLError, HTTPError

interpreter_path = sys.executable
pip_bootstrap_url = "https://bootstrap.pypa.io/get-pip.py"

dependencies = [
    "pyyaml"
]


def can_pip():
    try:
        import pip
    except ImportError, e:
        print(e)
        exit(-2)
    return True


def get_pip():
    try:
        f = urlopen(pip_bootstrap_url)
        with open(os.path.basename(pip_bootstrap_url), 'wb') as handle:
            handle.write(f.read())
    except HTTPError, e:
        print(e)
        exit(-3)
    except URLError, e:
        print(e)
        exit(-4)
    try:
        subprocess.check_call([interpreter_path, os.path.basename(pip_bootstrap_url)])
    except subprocess.CalledProcessError, e:
        print(e)
        exit(-5)


def install_dependencies():
    try:
        for dep in dependencies:
            subprocess.check_call(interpreter_path + " -m pip install %s --user" % dep)
    except subprocess.CalledProcessError, e:
        print(e)
        exit(-6)

if __name__ == '__main__':
    app = argparse.ArgumentParser()
    app.add_argument("-c", "--can-pip", action="store_true", default=False)
    app.add_argument("-g", "--get-pip", action="store_true", default=False)
    app.add_argument("-i", "--install_dependencies", action="store_true", default=False)

    args = app.parse_args()
    if args.can_pip:
        can_pip()
    if args.get_pip:
        get_pip()
    if args.install_dependencies:
        install_dependencies()
