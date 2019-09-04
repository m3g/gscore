#!/usr/bin/tclsh
#
# This script iterativelly checks if some process
# of some user is running. It continues to run until
# the process is not found anymore, in which case
# it quits. The checking is performed at each 0.1s, by default.
#
# It is used to wait for a process to finish to
# start a new one automatically.
#
# Run with:  wait_for.tcl user_name process_name [number_of_processes] [delay]
#
# Where user_name is the name of the owner of the process
# and the process_name is the process name, or PID, or
# anything that identifies univocally the process 
# of interest. number_of_processes is the maximum number
# of processes of that type allowed. That is, if the scripts
# finds at least number_of_processes processes running, then
# it will continue waiting. 
#
# Example:
#
# wait_for.tcl leandro namd2 4
#
# The script will run while there are at least 4 processes
# of user "leandro" that match "namd2". If one of the 4 processes
# finishes, the script exits. 
#
# L. Martinez, Sep 2, 2009.
#

# Just parse the command line 

proc parse_line { line char } {
  set line [ split $line $char ]
  set i 0
  set words(1) "#none"
  set words(2) "#none"
  set words(3) "#none"
  set words(4) "#none"
  foreach word $line { 
    if { $word > " " & $word != $char } { incr i; set words($i) $word }
  }
  array get words 
}

# Returns 1 if "n" processes that match the pattern are running, 0 otherwise
# Obs: Ignore the wait_for.tcl processes.

proc check_proc { user process n_processes } {
  catch { exec ps -eo user:50,cmd } ps  
  if { [ string first "wait_for" $ps ] != -1 } { 
    set ps [ split $ps "\n" ]
    set check_proc 0
    foreach line $ps {  
      if { [ string first $user $line ] != -1 &
           [ string first $process $line ] != -1 &
           [ string first "wait_for" $line ] == -1 } {
        incr check_proc 1     
      }
    }
  }
  if { $check_proc >= $n_processes } {
    return 1
  } else {
    return 0
  }
}

# The actual script

array set arg [ parse_line $argv " " ] 
set user $arg(1)
set process $arg(2)

if { $user == "#none" | $process == "#none" } {
  puts " Run with: wait_for.tcl user process_name \[n\]"
  exit
}

set n_processes 1
if { $arg(3) != "#none" } {  
  set n_processes $arg(3)
}

set time 0.1
if { $arg(4) != "#none" } {
  set time $arg(4)
}
set wait_check [ expr int(1000*$time) ]
set wait_confirm [ expr int($wait_check / 2) ]

set running 1
while { $running == 1 } {
  set running [ check_proc $user $process $n_processes ] 

# If the process is running wait 30 s

  if { $running } {
    after $wait_check

# If not, wait 15 seconds and check again (this is only
# because the process may be starting at the exact moment
# of the checking, just to be sure.

  } else {

    after $wait_confirm
    set running [ check_proc $user $process $n_processes ]
  
# Again, if the process is running wait 30 s and cycle

    if { $running } {  
      after $wait_check

# But if the process is not running, exit

    } else { 
      exit 
    }
  }
}


