#### SMILESBOX

##### Information:

A way of taking a SMILES string and turning it into a simple simulation cell 


##### Installation:

`pip install .`

##### Usage (CLI):

`smilesbox-generate -s "CC" --box --dimensions 10 10 10 --save`


##### Output:

```
SMILES:

 CC

FORMULA:

 C2H6

CELL:

[[10.  0.  0.]
 [ 0. 10.  0.]
 [ 0.  0. 20.]]

                                                  
                                                  
                                                  
                                                  
                 H                                
                                                  
                                                  
                  |          H                    
                  |                               
        H   ______|         /                     
                 _/\___   _/                      
                /      \_/______  H               
               /         |                        
             H           |                        
                         |                        
                                                  
                         H                        
                                        
```

![ethane molecule](https://github.com/badw/smilesbox/blob/devel/example/POSCAR.png?raw=true)


