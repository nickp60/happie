from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
record = next(SeqIO.parse(sys.argv[1], "genbank"))
print(record)
gd_diagram = GenomeDiagram.Diagram(record.id)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()
contig = SeqFeature(FeatureLocation(0, len(record.seq)))
gd_feature_set.add_feature(contig, sigil="ARROW",
                           color="black", label=True,
                           name="1",
                           arrowshaft_height=1.0,
                           label_size = 14, label_angle=0)
gd_feature_set.add_feature(contig, sigil="ARROW",
                           color="black", label=True,
                           name="1",
                           arrowshaft_height=1.0,
                           label_size = 14, label_angle=0)

for feature in record.features:
    if feature.type != "gene":
        #Exclude this feature
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue
    gd_feature_set.add_feature(feature, sigil="ARROW",
                               color=color, label=True,
                               label_size = 14, label_angle=0)

# #I want to include some strandless features, so for an example
# #will use EcoRI recognition sites etc.
# for site, name, color in [("GAATTC","EcoRI",colors.green),
#                           ("CCCGGG","SmaI",colors.orange),
#                           ("AAGCTT","HindIII",colors.red),
#                           ("GGATCC","BamHI",colors.purple)]:
#     index = 0
#     while True:
#         index  = record.seq.find(site, start=index)
#         if index == -1 : break
#         feature = SeqFeature(FeatureLocation(index, index+len(site)))
#         gd_feature_set.add_feature(feature, color=color, name=name,
#                                    label=True, label_size = 10,
#                                    label_color=color)
#         index += len(site)

gd_diagram.draw(format="linear", pagesize='A4', fragments=4,
                start=0, end=len(record))
gd_diagram.write("plasmid_linear_nice.pdf", "PDF")
gd_diagram.write("plasmid_linear_nice.eps", "EPS")
gd_diagram.write("plasmid_linear_nice.svg", "SVG")

gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                start=0, end=len(record), circle_core = 0.5)
gd_diagram.write("plasmid_circular_nice.pdf", "PDF")
gd_diagram.write("plasmid_circular_nice.eps", "EPS")
gd_diagram.write("plasmid_circular_nice.svg", "SVG")
