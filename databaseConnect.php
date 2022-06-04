<?php

// note must have ssh tunnel setup to forward port 8888 to the server:
// ssh -Nf -L 8888:localhost:3306 channagiri@knot.math.usf.edu
function tryConnectDatabase($database) {
  if ( $_SERVER['SERVER_NAME'] == 'localhost' ) {
    // running on the test computer
    // note must have a ssh tunnel setup on your computer to forward port 8888 to the server:
    // ssh -N -L 8888:localhost:3306 <username>@knot.math.usf.edu
    $link = mysqli_connect( 'localhost:8888', 'web', 'RZhRwsau6HZrMUXf', $database );
  } else {
    // running on the server
    $link = mysqli_connect( 'localhost', 'web', 'RZhRwsau6HZrMUXf', $database );
  }
  if ( mysqli_connect_errno() ) {
    printf( "Connect failed: %s\n", mysqli_connect_error() );
    exit();
  }
  return $link;  
}

?>