## Channel IDs used by HiggsSignals

These tables list the production and decay mode identifiers used in
HiggsSignals. A rate ID string is build out of one production mode ID `p` and one
decay mode ID `d` as `"p.d"`.

The production modes are:

 identifier | production mode
------------|-----------------
 0          | no production mode
 1          | single Higgs production, \f$ pp \to \phi \f$ (includes gluon fusion and bb-associated)
 2          | vector boson fusion, \f$ pp \to \phi \f$
 3          | W-boson associated Higgs production, \f$ pp \to W^\pm \phi \f$
 4          | Z-boson associated Higgs production, \f$ pp \to Z \phi \f$
 5          | tt-associated Higgs production, \f$ pp \to t\bar{t} \phi \f$
 6          | gluon fusion Higgs production, \f$ gg \to \phi \f$
 7          | bb-associated Higgs production, \f$ gg \to b\bar{b} \phi \f$/\f$ bb \to \phi \f$ (4FS/5FS)
 8          | single top associated Higgs production (t-channel), \f$ pp \to t\phi \f$
 9          | single top associated Higgs production (s-channel), \f$ pp \to t\phi \f$
 10         | quark-initiated Z-boson associated Higgs production, \f$ q\bar{q} \to Z\phi \f$
 11         | gluon-initiated Z-boson associated Higgs production, \f$ gg \to Z\phi \f$
 12         | single top and W-boson associated Higgs production, \f$ gb \to tW^\pm \phi \f$ (not yet implemented)

The decay modes are:

 identifier | decay mode
------------|-----------------
 0          | no decay mode
 1          | \f$ \phi \to \gamma\gamma \f$
 2          | \f$ \phi \to W^+W^- \f$
 3          | \f$ \phi \to ZZ \f$
 4          | \f$ \phi \to \tau^+\tau^- \f$
 5          | \f$ \phi \to b\bar{b} \f$
 6          | \f$ \phi \to Z\gamma \f$
 7          | \f$ \phi \to c\bar{c} \f$
 8          | \f$ \phi \to \mu^+\mu^- \f$
 9          | \f$ \phi \to gg \f$
 10         | \f$ \phi \to s\bar{s} \f$
 11         | \f$ \phi \to t\bar{t} \f$
