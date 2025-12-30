# generate_license.py
# Use this tool to produce a license.dat file for a given machine.
# Keep this script and your FERNET_KEY private.

import json
import uuid
import time
from cryptography.fernet import Fernet

# ====== CONFIGURE THIS VALUE (same key used in license_check.py) ======
FERNET_KEY = b'WbuZAYbiuCJQNLrQslxPv4a1WeLpna_R7MpDMej03o8='  # <--- same key as license_check.py
# ===================================================================
cipher = Fernet(FERNET_KEY)

def get_local_machine_id():
    return hex(uuid.getnode())

def make_license(name, email, machine_id, expiry="2026-12-31", allowed_runs=9999, license_key=None):
    data = {
        "name": name,
        "email": email,
        "license_key": license_key or "AUTO-" + str(int(time.time())),
        "machine_id": machine_id,
        "expiry": expiry,            # format YYYY-MM-DD
        "allowed_runs": int(allowed_runs)
    }
    enc = cipher.encrypt(json.dumps(data).encode("utf-8"))
    return enc

if __name__ == "__main__":
    print("Generating license.dat")
    name = input("Customer name: ").strip() or "Demo"
    email = input("Customer email: ").strip() or "demo@example.com"
    print("If you have a machine ID from the user, paste it now. Otherwise press Enter to use this machine's ID.")
    mid = input("Machine ID (hex): ").strip()
    if not mid:
        mid = get_local_machine_id()
        print(f"Using local machine id: {mid}")
    expiry = input("Expiry (YYYY-MM-DD) [default 2026-12-31]: ").strip() or "2026-12-31"
    runs = input("Allowed runs [default 9999]: ").strip() or "9999"
    boxed = make_license(name, email, mid, expiry=expiry, allowed_runs=int(runs))
    with open("license.dat", "wb") as fh:
        fh.write(boxed)
    print("license.dat written. Send this encrypted file to the user (do not share FERNET_KEY).")
