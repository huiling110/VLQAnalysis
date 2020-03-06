# VLQAnalysis: Analysis code to Search for Vector Like Quarks in an all Hadronic Final State

## Dependencies

This code requires the binary file of [BEST](https://gitlab.cern.ch/boostedeventshapetagger/BEST). This code is being
developed with `CMSSW_10_2_8`

## Repository Structure

BESTAnalyzer is an ntuplizer that makes root trees that output a BEST score.

## Some useful git commands

An example on how to add a remote and track a branch
```bash
git remote add reyer https://gitlab.cern.ch/rband/VLQAnalysis.git # this can be used to look at other peoples changes
git fetch -p --all #this gets all the branch information fro the remotes
git checkout -b ReyerMaster -t reyer/master #this makes a local branch 'ReyerMaster' that tracks Reyer's master branch
```

It is useful to make feature and bugfix branches whenever making changes.

Before pushing changes to be merged, it is useful to do a re-base.

```bash
# This rebases to Reyer master and adds your changes on top
git fetch -p --all
git checkout ReyerMaster
git pull
git checkout feature/myfeaturebranch
git rebase -i ReyerMaster
# then resolve any merge issues and continue with the rebase
git push
```

If after you pushed, someone else added more changes, then you will need to do another rebase and force push.
