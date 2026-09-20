"""screen/ — kill doomed candidates early. It does NOT rank survivors.

BUILD_PLAN §8.3 says what is reachable: "an open, honestly-labelled harness that
ranks candidates, kills doomed ones early, and shows its work." Those are two
jobs with very different evidence requirements, and conflating them is how a
screen starts producing numbers it has not earned.

    KILLING a candidate needs only to show it cannot work, or cannot be told
    apart from something that already exists, or breaks the model. Each of those
    is answerable today.

    RANKING survivors by predicted effect size needs the model to generalise to
    a mechanism it has never seen. `backtest/lomo.py` measures exactly that and
    reports 45.9pp against a 12.3pp predict-the-mean null. Until that inverts, a
    predicted effect size is not information.

So this package implements the first and REFUSES the second. `rank_candidates`
exists only to raise with the current LOMO number in the message.
"""

# Deliberately no eager re-export from `screen.kill_filter`. Importing it here
# makes `python -m screen.kill_filter` emit a RuntimeWarning about the module
# already being in sys.modules, and a warning nobody can act on trains people to
# ignore warnings. Import from the submodule:
#
#     from screen.kill_filter import screen, enumerate_candidates
