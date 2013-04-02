BEGIN{
	ORS="";
	ir=0;
	print "      double ";
	}
END{	
	print ";\n"
}	
	
	
$1~/[t]/ {
		ir++;
		if(ir>1)
		{
			if(ir<=11)
			{
				{print "," ;}
			}
		}	
		print $1
	        if(ir>10){print ";\n      double "; ir=0;}
			
         }
