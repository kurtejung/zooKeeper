# zooKeeper
RAA and RpA Plotting Utilities

The Zookeeper utilities allow users to plot all combinations of RpA and RAA plots from the CMS Collaboaration.  To use, simply download the git repository and modify (or add your own) the plot.txt file.
There can be either one column of keywords (for a single-panel plot) or two columns of keywords (for a two-panel plot), where the left column will be on the left plot and the right column will be on the right plot.

The standard structure of the plot.txt layout is:

ItemA <- tab -> ItemB </br>
ItemC </br>
ItemD </br>

The code allows you to put as many objects as you like on each panel of the plot.  If you wish to create a plot where column B is longer than column A, you'll need to literally write "placeholder" to properly load objects onto the right plot, i.e.

ItemA   < - tab ->    ItemB </br>
placeholder <- tab -> ItemC </br>
placeholder <- tab -> ItemD </br>

The current keywords are listed in <a href=https://github.com/kurtejung/zooKeeper/blob/master/drawZoo.C#L46>Line 46 of drawZoo.C</a>, however they are also listed here for convienience:
"HadronRAA", "HadronRpA", "InclJetRpA", "InclJetRAA", "BJetRAA", "PhotonRAA", "ZRAA", "WRAA", "BJpsiRAA", "BMesonRpA", "BJetRpA", "CJetRpA"

Finally, objects are drawn in the order of the input file, so the plots at the "back" of the figure should be first in the input file column, and vice-versa.

To add additional RpA or RAA data, please contact me
