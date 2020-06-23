#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/madman': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'AdapterRemoval': ['v_adapterremoval.txt', r"AdapterRemoval ver. (\S+)"],
    'fastp': ['v_fastp.txt', r"fastp (\S+)"],
    'MEGAHIT': ['v_megahit.txt', r"MEGAHIT v(\S+)"],
    'SPAdes': ['v_spades.txt', r"SPAdes genome assembler v(\S+)"],
    'QUAST': ['v_quast.txt', r"QUAST v(\S+)"],
    'Pydamage': ['v_pydamage.txt', r"pydamage, version (\S+)"],
    'DamageProfiler': ['v_damageprofiler.txt', r"DamageProfiler v(\S+)"],
    'prokka': ['v_prokka.txt', r"prokka (\S+)"]
}
results = OrderedDict()
results['nf-core/madman'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['AdapterRemoval'] = '<span style="color:#999999;\">N/A</span>'
results['fastp'] = '<span style="color:#999999;\">N/A</span>'
results['MEGAHIT'] = '<span style="color:#999999;\">N/A</span>'
results['SPAdes'] = '<span style="color:#999999;\">N/A</span>'
results['QUAST'] = '<span style="color:#999999;\">N/A</span>'
results['Pydamage'] = '<span style="color:#999999;\">N/A</span>'
results['DamageProfiler'] = '<span style="color:#999999;\">N/A</span>'
results['prokka'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Dump to YAML
print('''
id: 'software_versions'
section_name: 'nf-core/madman Software Versions'
section_href: 'https://github.com/nf-core/madman'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
