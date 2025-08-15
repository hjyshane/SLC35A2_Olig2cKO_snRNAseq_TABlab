gene_list <- list(Nsg2 = "Nsg2",
                  Aqp4 = "Aqp4",
                  Apoe = "Apoe",
                  P2ry12 = "P2ry12",
                  Prox1 = "Prox1",
                  Bcl11b = "Bcl11b",
                  Grik4 = "Grik4",
                  Neurod2 = "Neurod2",
                  Mbp = "Mbp",
                  Pdgfra = "Pdgfra",
                  Vip = "Vip",
                  Sst = "Sst",
                  Tbr1 = "Tbr1",
                  Satb2 = "Satb2",
                  Cux2 = "Cux2",
                  Rorb = "Rorb",
                  Slc17a7 = "Slc17a7",
                  Gad1 = "Gad1",
                  Vim = "Vim",
                  Ndnf = "Ndnf",
                  Reln = "Reln",
                  Dcn = "Dcn",
                  Dcx = "Dcx",
                  Tle4 = "Tle4",
                  Meis2 = "Meis2",
                  Lamp5 = "Lamp5",
                  Neurod6 = "Neurod6",
                  Meis2 = "Meis2",
                  Zfhx3 = "Zfhx3",
                  Lhx6 = "Lhx6",
                  # Meg3 = "Meg3",
                  Tubb2a = "Tubb2a"
)

genes_to_plot <- c(
    "Aqp4",  # Astrocytes
    "P2ry12",  # Microglia
    "Pdgfra",  # OPCs
    "Mbp",  # COPMFOL (Myelin-associated)
    "Dcn",  # CA1
    "Matn2",  # CA1
    "Grik4",  # CA3
    "Prox1",  # DG in general CA1-
    "Calb1",  # ExcitatoryNeuronsMatureDG
    "Dcx",  # ExcitatoryNeuronsImmatureDG
    "Satb2",  # cortical neuron
    "Cux2",  # L234ExcitatoryNeurons
    "Bcl11b",  # L5ExcitatoryNeurons CA1+
    "Rorb", # L4/5ExcitatoryNeurons
    "Etv1",  # L5ExcitatoryNeurons
    "Tbr1",  # L56ExcitatoryNeurons
    "Foxp2",  # L6
    "Gad2", # inhibitory neurons
    "Vip",  # CGE-derived cells
    "Meis2",  # LGE-derived interneurons
    "Sst",  # MGE-derived interneurons
    "Ndnf",  # NdnfRelnInhibitoryInterneurons
    "Lamp5",  # Lamp5 Positive
    "Col1a1",  # Vascular Cells
    "Vim"     # Marker for progenitor or stem-like cells
)

OligMarkers <- c(
    'Pdgfra', 'Cspg4', 'Olig1', 'Olig2', 'Sox10', 'Sox6',            # OPC markers
    'Ptprz1', 'Lhfpl3', 'Bmp4', 'Neu4', 'Id2', 'Id4', 'Tcf7l2',      # COP markers, Sox10 as well.
    'Enpp6', 'Mog', 'Plp1', 'Prom1', 'Tspan2', 'Ctps1', 'Grp17',     # NFOL markers
    'Mobp', 'Mag', 'Mal', 'Opalin', 'Serinc5', 'Cnp', 'Mbp', 'Myrf', # MFOL markers, Mog, Plp1 as well.
    'Cnp', 'Klk6', 'Apod', 'Trf', 'Pmp22'                            # Mature_OL markers, Mobp, Mbp, Mag, Opalin,  Plp1 as well
    )


OligMarkers_selected <- c(
    'Pdgfra', 'Olig2', 'Sox6',
    'Ptprz1', 'Neu4', 'Sox10',
    'Enpp6', 'Mog', 'Plp1', 
    'Mag', 'Mbp', 'Cnp'
)
