#!/bin/bash

bgzip OE19_DNMTi_HDACi_bedMethyl.bed
tabix -p bed OE19_DNMTi_HDACi_bedMethyl.bed.gz

