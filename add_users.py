from werkzeug.security import generate_password_hash
import shelve
import sys

username = sys.argv[1]
password = sys.argv[2]
action = sys.argv[3]

username_shelve = shelve.open("users.shelve", writeback=True)

if action == "add":
    username_shelve[username] = generate_password_hash(password)
elif action == "remove":
    del username_shelve[username]
