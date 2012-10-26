proc loadmultiple { filespec } { 
  foreach filename [glob $filespec] { 
    mol new $filename 
  } 
} 
