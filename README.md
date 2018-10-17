# FLASH5

## Git/Testing Workflow

The current rules for collaborating via git are as follows
1.  Base all feature branches off of the master branch.
2.  When development on a feature branch is finished, the feature branch is
first merged into the development branch.  The merge can be done directly by the
author or by a pull request assigned to a different developer.
3.  When the author has successfully executed on the development branch all 
manual tests related to changes made in the feature branch and believes that the
feature branch is ready for inclusion in FLASH5, a pull request is issued for
merging the feature branch into the staged branch.
    * One can trigger a code review at this point if so desired.
    * If all changes made in the feature branch are made to a portion of the code
that is not yet under automated testing and if it is believed that the changes can not break
any automated tests, then the pull request can be managed fully by the creator.
Else, the pull request should be assigned to Klaus or Jared.
4.  If the pull request is accepted, then the merge can be made and the person
that executes the merge must launch the compute001/Jenkins staged test suite ([Flash-Staged](https://jenkins-ci.mcs.anl.gov/job/Flash-Staged/)).
5.  If all tests in Flash-Staged have passed, then the person who performed the
merge into the staged branch issues a pull request from the feature branch to
the master branch.  This same individual accepts the request and subsequently
launches the compute001/Jenkins master test suite ([Flash-Master](https://jenkins-ci.mcs.anl.gov/job/Flash-Master/)).
6.  The person who triggered the Flash-Master test run follows up to confirm
that the master branch was not broken by the feature branch.
7.  To avoid accidentally including in the master branch commits that were
merged into the development or staged branches and subsequently abandoned, the
development and staged branches should never be merged into any other branch.
Note that when a merge conflict arises when merging a feature branch into one of
these infinite lifetime branches, attempts to resolve the conflicts using the
GitHub web interface will break this rule.
8.  Do not rebase a feature branch that has already been pushed to the GitHub
repository.

Note that neither Flash-Staged nor Flash-Master is being run on a regular basis
nor automatically in response to work done in GitHub.

In its present form, FlashTest does not allow for simultaneous execution of test
suites.  Therefore, when launching test suite runs on Jenkins, please confirm
first that no other runs are already in progress.  This can be done by consulting the [main Jenkins page](https://jenkins-ci.mcs.anl.gov), which shows all Jenkins runs in progress in the Build Executor Pane.

We presently run on a weekly basis the small subset of AMReX-only tests ([Flash-Staged-Intel](https://jenkins-ci.mcs.anl.gov/job/Flash-Staged-Intel/)) that are built in debug mode with the Intel compiler suite.

