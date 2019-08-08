import requests
import os
import json

payload = {
    'name': os.environ['name'],
    'institution': os.environ['institution'],
    'email': os.environ['email'],
    'product': 'htmd',
    'city': os.environ['city'],
    'country': os.environ['country']
}

product = 'htmd'

# try:
r = requests.post("https://www.acellera.com/licensing/htmd/register.php", params=payload)
ret = json.dumps(r.content.decode("ascii"))

print(ret)

prefix = os.path.join(os.path.expanduser('~'), '.htmd')
if not os.path.exists(prefix):
    os.mkdir(prefix)
prefix = os.path.join(prefix, '.registered-' + product)
if not os.path.exists(prefix):
    os.mkdir(prefix)
regfile = os.path.join(prefix, "registration")
fh = open(regfile, "w")
fh.write(r.content.decode("ascii"))
fh.close()
print("")
