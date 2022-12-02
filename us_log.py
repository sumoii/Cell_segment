import os,sys,time
import logging

logger = logging.getLogger(__name__)
"""
logger.setLevel(level = logging.INFO)
handler = logging.FileHandler("tmp.tmp"+"_log.txt")
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
"""

def log_setting(prefix):
    #logger = logging.getLogger(__name__)
    logger.setLevel(level = logging.INFO)
    handler = logging.FileHandler(prefix+"_log.txt")
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

