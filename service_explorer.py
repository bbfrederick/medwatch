#!/usr/bin/env python
	
"""
Service Explorer
----------------

An example showing how to access and print out the services, characteristics and
descriptors of a connected GATT server.

Created on 2019-03-25 by hbldh <henrik.blidh@nedomkull.com>

"""

import sys
import platform
import asyncio
import logging

from bleak import BleakClient

ADDRESS = (
    "24:71:89:cc:09:05"
    if platform.system() != "Darwin"
    else "EC155497-AF39-4C35-86FB-326F4226F45E"
)
if len(sys.argv) == 2:
    ADDRESS = sys.argv[1]


async def run(address, debug=False):
    log = logging.getLogger(__name__)
    if debug:
        import sys

        log.setLevel(logging.DEBUG)
        h = logging.StreamHandler(sys.stdout)
        h.setLevel(logging.DEBUG)
        log.addHandler(h)

    async with BleakClient(address) as client:
        log.info(f"Connected: {client.is_connected}")

        for service in client.services:
            log.info(f"[Service] {service}")
            for char in service.characteristics:
                print(char.properties)
                if "read" in char.properties:
                    try:
                        value = bytes(await client.read_gatt_char(char.uuid))
                        log.info(
                            f"\t[Characteristic] {char} ({','.join(char.properties)}), Value: {value}"
                        )
                    except Exception as e:
                        log.error(
                            f"\t[Characteristic] {char} ({','.join(char.properties)}), Value: {e}"
                        )

                else:
                    value = None
                    log.info(
                        f"\t[Characteristic] {char} ({','.join(char.properties)}), Value: {value}"
                    )

                for descriptor in char.descriptors:
                    try:
                        value = bytes(
                            await client.read_gatt_descriptor(descriptor.handle)
                        )
                        log.info(f"\t\t[Descriptor] {descriptor}) | Value: {value}")
                    except Exception as e:
                        log.error(f"\t\t[Descriptor] {descriptor}) | Value: {e}")


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    loop.set_debug(True)
    loop.run_until_complete(run(ADDRESS, True))
