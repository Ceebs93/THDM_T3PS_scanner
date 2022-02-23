from subprocess import run, PIPE

def keyFromLine(line):
    tableoffset = 83
    keylen = 35
    return line[tableoffset : tableoffset + keylen].strip()


def keySort(keylist):
    keys = list(set(filter(lambda x: len(x) > 0, keylist)))
    keys.sort(key=lambda x: x.split(":")[-1])
    return keys


res = run(["../AllAnalyses"], stdout=PIPE, universal_newlines=True)
with open("HBAnalyses.txt", "w") as allAnalyses:
    allAnalyses.write(res.stdout)
lines = res.stdout.split("\n")
lepKeys = keySort([keyFromLine(x) for x in lines if " LEP " in x])
tevKeys = keySort(
    [keyFromLine(x) for x in lines if " D0 " in x or " CDF " in x or " TCB " in x]
)
atlasKeys = keySort([keyFromLine(x) for x in lines if " ATL " in x])
cmsKeys = keySort([keyFromLine(x) for x in lines if " CMS " in x])


def getInspireBibtex(texkey):
    import requests

    inspire_query_params = {"q": "texkey=" + texkey, "format": "bibtex"}
    inspire_endpoint = "https://inspirehep.net/api/literature"
    res = requests.get(inspire_endpoint, params=inspire_query_params)
    bibtex = res.text
    return bibtex


with open("HBAnalyses.bib", "w") as bibfile:
    for texkey in lepKeys + tevKeys + atlasKeys + cmsKeys:
        bibfile.write(
            getInspireBibtex(texkey)
            .replace("\u2212", "-")
            .replace(" &", r" \&")
            .replace(r"\\tau", r"\tau")
            .replace("find_paper", r"find\_paper")
        )
    bibfile.write(
        """
@inproceedings{CDFNotes,
    title = "unpublished CDF Notes",
    reportNumber  = "{CDF-9999, CDF-10010, CDF-10439, CDF-10500, CDF-10573, CDF-10574}",
    author        = "{CDF}",
    collaboration = "{CDF}",
}

@inproceedings{D0Notes,
    title = "unpublished D0 Notes",
    reportNumber  = "{D0-5845, D0-6183, D0-6305}",
    author        = "{D0}",
    collaboration = "{D0}",
}
"""
    )

texstring = (
    "HiggsBounds currently incorporates results from\n"
    "LEP~\\cite{{{lepcite}}},\n"
    "the Tevatron~\\cite{{{tevcite}, CDFNotes, D0Notes}},\n"
    "and the ATLAS~\\cite{{{atlascite}}}\n"
    "and CMS~\\cite{{{cmscite}}}\n"
    "experiments at the LHC\@.\n"
)

with open("HBAnalyses.tex", "w") as texfile:
    texfile.write(
        r"""\documentclass{scrartcl}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{hyperref}
\usepackage{textgreek}
\usepackage{amsmath}
\usepackage{multicol}
\usepackage[citestyle=numeric-comp,bibstyle=HBnumeric,sorting=none]{biblatex}
\addbibresource{HBAnalyses.bib}
\begin{document}
\section*{References for the analyses currently implemented in \texttt{HiggsBounds}}
"""
    )
    texfile.write(
        texstring.format(
            lepcite=", ".join(lepKeys),
            tevcite=", ".join(tevKeys),
            atlascite=", ".join(atlasKeys),
            cmscite=", ".join(cmsKeys),
        )
    )
    texfile.write(
        r"""
\begin{multicols}{2}[\printbibheading]
	\renewcommand*{\bibfont}{\small}
	\printbibliography[heading=none]
\end{multicols}
\end{document}"""
    )
