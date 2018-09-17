# FLASH5

[![Build status from Travis CI](https://img.shields.io/travis/com/mbdevpl/FLASH5/travis.svg)](https://travis-ci.com/mbdevpl/FLASH5)

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

## Tests on Master Branch

The following are the setup commands of the tests that are currently being run on the Master branch to confirm correctness of the functionality we consider functional in FLASH5.

#### Unit Tests
* unitTest/Eos/Helmholtz -auto +amrex -3d +noio
* unitTest/Eos/Helmholtz -auto +pm4dev -3d +noio
* unitTest/Gravity/Poisson3 -auto -2d +cylindrical +newmpole -debug -maxblocks=600 +noio +pm4dev -parfile=flash_2dcyl.par
* unitTest/Gravity/Poisson3 -auto -3d +newmpole +uhd -debug -maxblocks=550 -nxb=8 -nyb=8 -nzb=8 -gridinterpolation=monotonic
* unitTest/Grid/Amrex/TestFluxCorrection -auto -2d -nxb=8 -nyb=8 +noio +amrex
* unitTest/Grid/Amrex/TestFluxCorrection2 -auto -2d -nxb=8 -nyb=8 +noio +amrex
* unitTest/Grid/Amrex/TestInit -auto -2d -nxb=8 -nyb=4 +noio +amrex
* unitTest/Grid/Amrex/TestRefine -auto -2d -nxb=8 -nyb=8 +noio +amrex
* unitTest/Grid/Amrex/TestCyl2d -auto -2d -nxb=8 -nyb=4 +noio +amrex
* unitTest/Multigrid_Amrex -auto -3d +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM -maxblocks=1000

#### Regression Tests with Verified Baseline
* Sod -auto -2d -debug +uhd +ug +nofbs -parfile=test_pseudoug_2d.par
* Sod -auto -2d -debug -nxb=8 -nyb=8 +uhd +pm4dev -gridinterpolation=monotonic -parfile=test_amr_unsplit_2d.par
* Sod -auto -2d -debug -nxb=8 -nyb=8 -unit=physics/Hydro/HydroMain/simpleUnsplit/HLL +pm4dev -gridinterpolation=native -parfile=test_amr_2d.par
* Sod -auto -2d -debug -nxb=8 -nyb=8 +uhd +amrex -gridinterpolation=native +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=	test_amr_unsplit_2d.par
* Sod -auto -2d -debug -nxb=8 -nyb=8 -unit=physics/Hydro/HydroMain/simpleUnsplit/HLL +amrex -gridinterpolation=native +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test_amr_2d.par
* Sedov -auto -2d -debug +uhd +ug +nofbs -parfile=	test_pseudoug_2d.par
* Sedov -auto -3d -debug -nxb=8 -nyb=8 -nzb=8 +uhd +pm4dev -gridinterpolation=native -parfile=test_amr_unsplit_3d.par
* Sedov -auto -3d -debug -nxb=8 -nyb=8 -nzb=8 +uhd +amrex -gridinterpolation=native +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test_amr_unsplit_3d.par
* Sedov -auto -2d +cylindrical -debug -nxb=16 -nyb=16 +uhd +pm4dev -gridinterpolation=monotonic DoAnalytical=True -parfile=	test_amr_cyl_2d.par
* Sedov -auto -2d +cylindrical -debug -nxb=16 -nyb=16 +uhd +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM DoAnalytical=True -parfile=test_amr_cyl_2d.par
* Sedov -auto -2d -debug -nxb=8 -nyb=8 +uhd +pm4dev -gridinterpolation=native -parfile=test_amr_unsplit_2d.par
* Sedov -auto -2d -debug -nxb=8 -nyb=8 -unit=physics/Hydro/HydroMain/simpleUnsplit/HLL +pm4dev -gridinterpolation=native -parfile=test_amr_2d.par
* Sedov -auto -2d -debug -nxb=8 -nyb=8 +uhd +amrex -gridinterpolation=native +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test_amr_unsplit_2d.par
* Sedov -auto -2d -debug -nxb=8 -nyb=8 -unit=physics/Hydro/HydroMain/simpleUnsplit/HLL +amrex -gridinterpolation=native +serialio -unit=IO/IOMain/hdf5/serial/AM -parfile=test_amr_2d.par
* Cellular -auto -2d -debug +a13 +uhd +pm4dev -gridinterpolation=monotonic -parfile=test_amr_2d.par

#### Regression Tests with Unverified Baseline
* DustCollapse -auto -3d +cartesian +Mode1 +serialIO +uhd +newMpole -debug -parfile=test_3dcar.par
* DustCollapse -auto -3d +cartesian +Mode3 +serialIO +uhd +newMpole -debug -parfile=test_3dcar.par
* DustCollapse -auto -2d +cylindrical +Mode1 +serialIO +uhd +newMpole -debug -parfile=	test_2dcyl.par
* DustCollapse -auto -2d +cylindrical +Mode3 +serialIO +uhd +newMpole -debug -parfile=test_2dcyl.debug.par
* DustCollapse -auto -1d +spherical +Mode1 +serialIO +uhd +newMpole -debug -parfile=test_1dsph.par
* DustCollapse -auto -1d +spherical +Mode3 +serialIO +uhd +newMpole -debug -parfile=test_1dsph.debug.par
* IsentropicVortex -auto -2d -debug +uhd +amrex +serialIO -unit=IO/IOMain/hdf5/serial/AM withParticles=TRUE
