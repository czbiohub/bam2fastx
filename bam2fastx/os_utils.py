import os

# Set input language USA unicode encoding setting
# Necessary because click assumes ascii input unless otherwise specified
# https://click.palletsprojects.com/en/7.x/python3/
unicode_usa = 'en_US.utf-8'
os.environ['LC_LANG'] = unicode_usa
os.environ["LC_ALL"] = unicode_usa
