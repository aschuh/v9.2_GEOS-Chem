	  o  F   k820309    V          13.1        jOaT                                                                                                           
       bpch2_mod.F BPCH2_MOD              OPEN_BPCH2_FOR_READ OPEN_BPCH2_FOR_WRITE BPCH2_HDR BPCH2 READ_BPCH2 GET_MODELNAME GET_NAME_EXT GET_NAME_EXT_2D GET_RES_EXT GET_HALFPOLAR gen@GET_TAU0                      @                              
       FINDFREELUN I_AM_UNOPENED                                                        u #GET_TAU0_6A    %         @    @X                                               
       #GET_TAU0_6A%DBLE    #GET_TAU0_6A%PRESENT    #MONTH    #DAY    #YEAR    #HOUR    #MIN 	   #SEC 
                 @                                 DBLE               @                                 PRESENT           
  @@                                                   
  @@                                                   
  @@                                                   
 @@                                                   
 @@                              	                     
 @@                              
           %         @      !                                                    #FINDFREELUN%TRIM    #B                  @                                 TRIM           
  @                                         %         @      !                                                      #N                @                                          #         @                                                     #OPEN_BPCH2_FOR_READ%PRESENT    #OPEN_BPCH2_FOR_READ%TRIM    #IUNIT    #FILENAME    #TITLE                  @                                 PRESENT               @                                 TRIM           
  @@                                                   
  @@                                                 1           F @@                                  P                       #         @                                                      #OPEN_BPCH2_FOR_WRITE%PRESENT    #OPEN_BPCH2_FOR_WRITE%TRIM    #IUNIT    #FILENAME    #TITLE                  @                                 PRESENT               @                                 TRIM           
  @@                                                   
  @@                                                 1            @@                                  P                       #         @                                                      #IUNIT    #TITLE              
  @@                                                   
   @                                  P                      #         @                                                       #IUNIT     #MODELNAME !   #LONRES "   #LATRES #   #HALFPOLAR $   #CENTER180 %   #CATEGORY &   #NTRACER '   #UNIT (   #TAU0 )   #TAU1 *   #RESERVED +   #NI ,   #NJ -   #NL .   #IFIRST /   #JFIRST 0   #LFIRST 1   #ARRAY 2             
  @@                                                    
   @                             !                                     
   @                             "     	                
   @                             #     	                
   @                              $                     
   @                              %                     
   @                             &     (                                
   @                              '                     
   @                             (     (                                
   @                             )     
                
   @                             *     
                
   @                             +     (                                
   @                              ,                     
   @                              -                     
   @                              .                     
   @                              /                     
   @                              0                     
   @                              1                    
   @                             2                    	        p        5 � p        r -   p        5 � p        r ,   p          5 � p        r ,     5 � p        r -     5 � p        r .       5 � p        r ,     5 � p        r -     5 � p        r .                     #         @                                   3                	   #READ_BPCH2%PRESENT 4   #READ_BPCH2%TRIM 5   #FILENAME 6   #CATEGORY_IN 7   #TRACER_IN 8   #TAU0_IN 9   #IX :   #JX ;   #LX <   #ARRAY =   #QUIET >                 @                            4     PRESENT               @                            5     TRIM           
  @@                             6                    1           
  @@                             7                    1           
   @                              8                     
   @                             9     
                
   @                              :                     
   @                              ;                     
   @                              <                    D  @                             =                    	         p        5 � p        r ;   p        5 � p        r :   p          5 � p        r :     5 � p        r ;     5 � p        r <       5 � p        r :     5 � p        r ;     5 � p        r <                               
 @@                              >           $         @                                 ?                                    $         @                                @                                    $         @                                 A                                    $         @                                 B                                    %         @                                 C                               �         fn#fn    �   �   b   uapp(BPCH2_MOD    d  Z   J  INQUIREMOD    �  Q       gen@GET_TAU0      �      GET_TAU0_6A !   �  =      GET_TAU0_6A%DBLE $     @      GET_TAU0_6A%PRESENT "   E  @   a   GET_TAU0_6A%MONTH     �  @   a   GET_TAU0_6A%DAY !   �  @   a   GET_TAU0_6A%YEAR !     @   a   GET_TAU0_6A%HOUR     E  @   a   GET_TAU0_6A%MIN     �  @   a   GET_TAU0_6A%SEC '   �  m       FINDFREELUN+INQUIREMOD 1   2  =      FINDFREELUN%TRIM+INQUIREMOD=TRIM )   o  @   e   FINDFREELUN%B+INQUIREMOD )   �  W       I_AM_UNOPENED+INQUIREMOD +     @   e   I_AM_UNOPENED%N+INQUIREMOD $   F  �       OPEN_BPCH2_FOR_READ ,   �  @      OPEN_BPCH2_FOR_READ%PRESENT )   1  =      OPEN_BPCH2_FOR_READ%TRIM *   n  @   a   OPEN_BPCH2_FOR_READ%IUNIT -   �  L   a   OPEN_BPCH2_FOR_READ%FILENAME *   �  P   a   OPEN_BPCH2_FOR_READ%TITLE %   J  �       OPEN_BPCH2_FOR_WRITE -   �  @      OPEN_BPCH2_FOR_WRITE%PRESENT *   7	  =      OPEN_BPCH2_FOR_WRITE%TRIM +   t	  @   a   OPEN_BPCH2_FOR_WRITE%IUNIT .   �	  L   a   OPEN_BPCH2_FOR_WRITE%FILENAME +    
  P   a   OPEN_BPCH2_FOR_WRITE%TITLE    P
  ^       BPCH2_HDR     �
  @   a   BPCH2_HDR%IUNIT     �
  P   a   BPCH2_HDR%TITLE    >  &      BPCH2    d  @   a   BPCH2%IUNIT     �  P   a   BPCH2%MODELNAME    �  @   a   BPCH2%LONRES    4  @   a   BPCH2%LATRES     t  @   a   BPCH2%HALFPOLAR     �  @   a   BPCH2%CENTER180    �  P   a   BPCH2%CATEGORY    D  @   a   BPCH2%NTRACER    �  P   a   BPCH2%UNIT    �  @   a   BPCH2%TAU0      @   a   BPCH2%TAU1    T  P   a   BPCH2%RESERVED    �  @   a   BPCH2%NI    �  @   a   BPCH2%NJ    $  @   a   BPCH2%NL    d  @   a   BPCH2%IFIRST    �  @   a   BPCH2%JFIRST    �  @   a   BPCH2%LFIRST    $  �  a   BPCH2%ARRAY    �  �       READ_BPCH2 #   �  @      READ_BPCH2%PRESENT     �  =      READ_BPCH2%TRIM $     L   a   READ_BPCH2%FILENAME '   _  L   a   READ_BPCH2%CATEGORY_IN %   �  @   a   READ_BPCH2%TRACER_IN #   �  @   a   READ_BPCH2%TAU0_IN    +  @   a   READ_BPCH2%IX    k  @   a   READ_BPCH2%JX    �  @   a   READ_BPCH2%LX !   �  �  a   READ_BPCH2%ARRAY !     @   a   READ_BPCH2%QUIET    �  X       GET_MODELNAME      X       GET_NAME_EXT     o  X       GET_NAME_EXT_2D    �  X       GET_RES_EXT      P       GET_HALFPOLAR 