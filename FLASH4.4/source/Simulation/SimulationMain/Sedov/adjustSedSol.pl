#! /usr/bin/perl -p
$E=1.0, $rho_up=1.0, $t=0.01;
$minRho_Flash = 0.001;
$minRho_Flash_sq = $minRho_Flash * $minRho_Flash;
$minRho_Flash_4th = $minRho_Flash_sq * $minRho_Flash_sq;
$S =  ($rho_up / ($E*$t*$t) )**0.2 ;
#s/^# SOLUTION:/# ADJUSTED SOLUTION: E=$E, rho_up=$rho_up, t=$t  =>  S = $S/ ;
s/^# SOLUTION:/# SOLUTION WITH DENSITY FLOOR:/ ;
unless (m/^#/) {
    s/^((?:\s+[\d.Ee+-]+){2}\s+)([\d.Ee+-]+)\s(?=\s)/  # modify the 3rd number, leaving others as they are
               "$1" . 
         sprintf("%14.9E",
                 ( ($2 < 10.0 * $minRho_Flash)  ?
                     sqrt(sqrt($2**4 + $minRho_Flash_4th)) : $2 )
                )
     /e
}
