	  !  M   k820309    V          13.1        iOaT                                                                                                           
       error_mod.F ERROR_MOD              ALLOC_ERR DEBUG_MSG ERROR_STOP GEOS_CHEM_STOP IS_SAFE_DIV IS_SAFE_EXP SAFE_DIV SAFE_EXP SAFE_LOG SAFE_LOG10 gen@IT_IS_NAN gen@IT_IS_FINITE gen@CHECK_VALUE                                                        u #NAN_FLOAT    #NAN_DBLE    %         @    @X                                                     #NAN_FLOAT%ISNAN    #VALUE                  @                                 ISNAN           
  @@                                  	      %         @    @X                                                     #NAN_DBLE%ISNAN    #VALUE                  @                                 ISNAN           
  @@                                  
                                                             u #FINITE_FLOAT    #FINITE_DBLE 
   %         @    @X                                                     #FINITE_FLOAT%FP_CLASS    #VALUE 	                 @                                 FP_CLASS           
  @@                             	     	      %         @    @X                           
                          #FINITE_DBLE%FP_CLASS    #VALUE                  @                                 FP_CLASS           
  @@                                  
                                                             u #CHECK_REAL_VALUE    #CHECK_DBLE_VALUE    #         @     @X                                                #CHECK_REAL_VALUE%TRIM    #CHECK_REAL_VALUE%REPEAT    #VALUE    #LOCATION    #VARNAME    #MESSAGE                  @                                 TRIM               @                                 REPEAT           
@ @@                                  	                
   @                                                     p          p            p                                    
  @@                                  �                                
  @@                                  �                      #         @     @X                                                #CHECK_DBLE_VALUE%TRIM    #CHECK_DBLE_VALUE%REPEAT    #VALUE    #LOCATION    #VARNAME    #MESSAGE                  @                                 TRIM               @                                 REPEAT           
@ @@                                  
                
   @                                                     p          p            p                                    
  @@                                  �                                
  @@                                  �                      #         @                                                      #ALLOC_ERR%PRESENT    #ALLOC_ERR%TRIM    #ARRAYNAME    #AS                  @                                 PRESENT               @                                 TRIM           
  @@                                                 1           
B @@                                         #         @                                                        #MESSAGE !             
   @                             !                    1 #         @                                  "                   #ERROR_STOP%TRIM #   #ERROR_STOP%REPEAT $   #MESSAGE %   #LOCATION &                 @                            #     TRIM               @                            $     REPEAT           
  @@                             %                    1           
  @@                             &                    1 #         @                                  '                                              %         @                                 (                          #IS_SAFE_DIV%MINEXPONENT )   #IS_SAFE_DIV%MAXEXPONENT *   #IS_SAFE_DIV%EXPONENT +   #IS_SAFE_DIV%PRESENT ,   #N -   #D .   #R4 /                 @                            )     MINEXPONENT               @                            *     MAXEXPONENT               @                            +     EXPONENT               @                            ,     PRESENT           
  @@                             -     
                
  @@                             .     
                
 @@                              /           %         @                                0                          #IS_SAFE_EXP%ABS 1   #X 2                 @                            1     ABS           
  @@                             2     
      %         @                                 3                   
       #SAFE_DIV%MINEXPONENT 4   #SAFE_DIV%MAXEXPONENT 5   #SAFE_DIV%EXPONENT 6   #SAFE_DIV%PRESENT 7   #N 8   #D 9   #ALT_NAN :   #ALT_OVER ;   #ALT_UNDER <                 @                            4     MINEXPONENT               @                            5     MAXEXPONENT               @                            6     EXPONENT               @                            7     PRESENT           
  @@                             8     
                
  @@                             9     
                
   @                             :     
                
 @@                             ;     
                
 @@                             <     
      %         @                                 =                   
       #SAFE_EXP%EXP >   #X ?   #ALT @                 @                            >     EXP           
  @@                             ?     
                
   @                             @     
      %         @                                 A                   
       #SAFE_LOG%LOG B   #X C   #ALT D                 @                            B     LOG           
  @@                             C     
                
   @                             D     
      %         @                                 E                   
       #SAFE_LOG10%LOG10 F   #X G   #ALT H                 @                            F     LOG10           
  @@                             G     
                
   @                             H     
         �         fn#fn    �   �   b   uapp(ERROR_MOD    i  ]       gen@IT_IS_NAN    �  p      NAN_FLOAT     6  >      NAN_FLOAT%ISNAN     t  @   a   NAN_FLOAT%VALUE    �  o      NAN_DBLE    #  >      NAN_DBLE%ISNAN    a  @   a   NAN_DBLE%VALUE !   �  c       gen@IT_IS_FINITE      v      FINITE_FLOAT &   z  A      FINITE_FLOAT%FP_CLASS #   �  @   a   FINITE_FLOAT%VALUE    �  u      FINITE_DBLE %   p  A      FINITE_DBLE%FP_CLASS "   �  @   a   FINITE_DBLE%VALUE     �  l       gen@CHECK_VALUE !   ]  �      CHECK_REAL_VALUE &     =      CHECK_REAL_VALUE%TRIM (   M  ?      CHECK_REAL_VALUE%REPEAT '   �  @   a   CHECK_REAL_VALUE%VALUE *   �  �   a   CHECK_REAL_VALUE%LOCATION )   `  P   a   CHECK_REAL_VALUE%VARNAME )   �  P   a   CHECK_REAL_VALUE%MESSAGE !    	  �      CHECK_DBLE_VALUE &   �	  =      CHECK_DBLE_VALUE%TRIM (   �	  ?      CHECK_DBLE_VALUE%REPEAT '   /
  @   a   CHECK_DBLE_VALUE%VALUE *   o
  �   a   CHECK_DBLE_VALUE%LOCATION )     P   a   CHECK_DBLE_VALUE%VARNAME )   S  P   a   CHECK_DBLE_VALUE%MESSAGE    �  �       ALLOC_ERR "   -  @      ALLOC_ERR%PRESENT    m  =      ALLOC_ERR%TRIM $   �  L   a   ALLOC_ERR%ARRAYNAME    �  @   a   ALLOC_ERR%AS    6  U       DEBUG_MSG "   �  L   a   DEBUG_MSG%MESSAGE    �  �       ERROR_STOP     f  =      ERROR_STOP%TRIM "   �  ?      ERROR_STOP%REPEAT #   �  L   a   ERROR_STOP%MESSAGE $   .  L   a   ERROR_STOP%LOCATION    z  a       GEOS_CHEM_STOP    �  �       IS_SAFE_DIV (   �  D      IS_SAFE_DIV%MINEXPONENT (   �  D      IS_SAFE_DIV%MAXEXPONENT %   6  A      IS_SAFE_DIV%EXPONENT $   w  @      IS_SAFE_DIV%PRESENT    �  @   a   IS_SAFE_DIV%N    �  @   a   IS_SAFE_DIV%D    7  @   a   IS_SAFE_DIV%R4    w  l       IS_SAFE_EXP     �  <      IS_SAFE_EXP%ABS      @   a   IS_SAFE_EXP%X    _  �       SAFE_DIV %   H  D      SAFE_DIV%MINEXPONENT %   �  D      SAFE_DIV%MAXEXPONENT "   �  A      SAFE_DIV%EXPONENT !     @      SAFE_DIV%PRESENT    Q  @   a   SAFE_DIV%N    �  @   a   SAFE_DIV%D !   �  @   a   SAFE_DIV%ALT_NAN "     @   a   SAFE_DIV%ALT_OVER #   Q  @   a   SAFE_DIV%ALT_UNDER    �  r       SAFE_EXP      <      SAFE_EXP%EXP    ?  @   a   SAFE_EXP%X      @   a   SAFE_EXP%ALT    �  r       SAFE_LOG    1  <      SAFE_LOG%LOG    m  @   a   SAFE_LOG%X    �  @   a   SAFE_LOG%ALT    �  v       SAFE_LOG10 !   c  >      SAFE_LOG10%LOG10    �  @   a   SAFE_LOG10%X    �  @   a   SAFE_LOG10%ALT 