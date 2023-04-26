import os
from ipyparallel import Client
      
ufile = os.getenv('IPYTHONDIR') + "/profile_" \
        + os.getenv('IPYTHON_PROFILE') \
        + "/security/ipcontroller-client.json"

print(ufile)

rc = Client(url_file=ufile)

rc.shutdown(hub=True)

