	  �&  a   k820309    V          13.1        jOaT                                                                                                           
       grid_mod.F90 GRID_MOD              CLEANUP_GRID COMPUTE_GRID GET_AREA_M2 GET_AREA_CM2 GET_BOUNDING_BOX GET_XEDGE GET_XMID GET_YEDGE GET_YEDGE_R GET_YMID GET_YMID_R GET_YMID_R_W GET_YSIN GET_XOFFSET GET_YOFFSET INIT_GRID ITS_A_NESTED_GRID SET_XOFFSET SET_YOFFSET                                                    
                            @                             
       #         @                                                     #ERROR_STOP%REPEAT    #ERROR_STOP%TRIM    #MESSAGE    #LOCATION                  @                                 REPEAT               @                                 TRIM           
   @                                                 1           
   @                                                 1 #         @                                                     #ALLOC_ERR%TRIM 	   #ALLOC_ERR%PRESENT 
   #ARRAYNAME    #AS                  @                            	     TRIM               @                            
     PRESENT           
   @                                                 1           
  @                                         #         @                                                       #CLEANUP_GRID%ALLOCATED                  @                                 ALLOCATED #         @                                                      #COMPUTE_GRID%SIN    #AM_I_ROOT    #I1    #I2    #J1    #J2    #JSP    #JNP    #L1    #L2    #DLON    #DLAT    #I_LO    #J_LO    #RC                                                 @                                 SIN           
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                  
   @                                                 
        p          5 � p        r    5 � p        r    n                                           1p          5 � p        r    5 � p        r    n                                          1p            5 � p        r    5 � p        r    n                                      1    5 � p        r    5 � p        r    n                                          1    5 � p 	       r    5 � p        r    n                                          1      5 � p        r    5 � p        r    n                                      1    5 � p        r    5 � p        r    n                                          1    5 � p 	       r    5 � p        r    n                                          1                                          
   @                                                 
        p          5 � p        r    5 � p        r    n                                           1p          5 � p        r    5 � p        r    n                                          1p            5 � p        r    5 � p        r    n                                      1    5 � p        r    5 � p        r    n                                          1    5 � p 	       r    5 � p        r    n                                          1      5 � p        r    5 � p        r    n                                      1    5 � p        r    5 � p        r    n                                          1    5 � p 	       r    5 � p        r    n                                          1                                           
   @                                                   
   @                                                   D  @                                          %         @                                                     
       #I     #J !   #L "             
   @                                                    
   @                              !                     
   @                              "           %         @                                 #                    
       #I $   #J %   #L &             
   @                              $                     
   @                              %                     
   @                              &           #         @                                   '                    #I1 (   #I2 )   #J1 *   #J2 +   #L ,   #COORDS -   #INDICES .             
   @                              (                     
   @                              )                     
   @                              *                     
   @                              +                     
   @                              ,                     
   @                             -                   
    p          p            p                                    D  @                              .                        p          p            p                          %         @                                 /                    
       #I 0   #J 1   #L 2             
   @                              0                     
   @                              1                     
   @                              2           %         @                                 3                    
       #I 4   #J 5   #L 6             
   @                              4                     
   @                              5                     
   @                              6           %         @                                 7                    
       #I 8   #J 9   #L :             
   @                              8                     
   @                              9                     
   @                              :           %         @                                 ;                    
       #I <   #J =   #L >             
   @                              <                     
   @                              =                     
   @                              >           %         @                                 ?                    
       #I @   #J A   #L B             
   @                              @                     
   @                              A                     
   @                              B           %         @                                 C                    
       #I D   #J E   #L F             
   @                              D                     
   @                              E                     
   @                              F           %         @                                 G                    
       #I H   #J I   #L J             
   @                              H                     
   @                              I                     
   @                              J           %         @                                 K                    
       #I L   #J M   #L N             
   @                              L                     
   @                              M                     
   @                              N           %         @                                 O                          #GET_XOFFSET%PRESENT P   #GLOBAL Q                 @                            P     PRESENT           
 @@                              Q           %         @                                 R                          #GET_YOFFSET%PRESENT S   #GLOBAL T                 @                            S     PRESENT           
 @@                              T           #         @                                   U                    #AM_I_ROOT V   #IM W   #JM X   #LM Y   #RC Z                                         
   @                              V                     
   @                              W                     
   @                              X                     
   @                              Y                     D  @                              Z            %         @                                 [                            #         @                                   \                    #X_OFFSET ]             
   @                              ]           #         @                                   ^                    #Y_OFFSET _             
   @                              _              �         fn#fn    �   �   b   uapp(GRID_MOD    �  @   J  CMN_GCTM_MOD    �  @   J  ERROR_MOD %   1  �       ERROR_STOP+ERROR_MOD 3   �  ?      ERROR_STOP%REPEAT+ERROR_MOD=REPEAT /   �  =      ERROR_STOP%TRIM+ERROR_MOD=TRIM -   <  L   e   ERROR_STOP%MESSAGE+ERROR_MOD .   �  L   e   ERROR_STOP%LOCATION+ERROR_MOD $   �  �       ALLOC_ERR+ERROR_MOD .   ^  =      ALLOC_ERR%TRIM+ERROR_MOD=TRIM 4   �  @      ALLOC_ERR%PRESENT+ERROR_MOD=PRESENT .   �  L   e   ALLOC_ERR%ARRAYNAME+ERROR_MOD '   '  @   e   ALLOC_ERR%AS+ERROR_MOD    g  d       CLEANUP_GRID '   �  B      CLEANUP_GRID%ALLOCATED      �       COMPUTE_GRID !     <      COMPUTE_GRID%SIN '   G  @   a   COMPUTE_GRID%AM_I_ROOT     �  @   a   COMPUTE_GRID%I1     �  @   a   COMPUTE_GRID%I2       @   a   COMPUTE_GRID%J1     G  @   a   COMPUTE_GRID%J2 !   �  @   a   COMPUTE_GRID%JSP !   �  @   a   COMPUTE_GRID%JNP     	  @   a   COMPUTE_GRID%L1     G	  @   a   COMPUTE_GRID%L2 "   �	  \  a   COMPUTE_GRID%DLON "   �  \  a   COMPUTE_GRID%DLAT "   ?  @   a   COMPUTE_GRID%I_LO "     @   a   COMPUTE_GRID%J_LO     �  @   a   COMPUTE_GRID%RC    �  e       GET_AREA_M2    d  @   a   GET_AREA_M2%I    �  @   a   GET_AREA_M2%J    �  @   a   GET_AREA_M2%L    $  e       GET_AREA_CM2    �  @   a   GET_AREA_CM2%I    �  @   a   GET_AREA_CM2%J    	  @   a   GET_AREA_CM2%L !   I  �       GET_BOUNDING_BOX $   �  @   a   GET_BOUNDING_BOX%I1 $     @   a   GET_BOUNDING_BOX%I2 $   Q  @   a   GET_BOUNDING_BOX%J1 $   �  @   a   GET_BOUNDING_BOX%J2 #   �  @   a   GET_BOUNDING_BOX%L (     �   a   GET_BOUNDING_BOX%COORDS )   �  �   a   GET_BOUNDING_BOX%INDICES    9  e       GET_XEDGE    �  @   a   GET_XEDGE%I    �  @   a   GET_XEDGE%J      @   a   GET_XEDGE%L    ^  e       GET_XMID    �  @   a   GET_XMID%I      @   a   GET_XMID%J    C  @   a   GET_XMID%L    �  e       GET_YEDGE    �  @   a   GET_YEDGE%I    (  @   a   GET_YEDGE%J    h  @   a   GET_YEDGE%L    �  e       GET_YEDGE_R      @   a   GET_YEDGE_R%I    M  @   a   GET_YEDGE_R%J    �  @   a   GET_YEDGE_R%L    �  e       GET_YMID    2  @   a   GET_YMID%I    r  @   a   GET_YMID%J    �  @   a   GET_YMID%L    �  e       GET_YMID_R    W  @   a   GET_YMID_R%I    �  @   a   GET_YMID_R%J    �  @   a   GET_YMID_R%L      e       GET_YMID_R_W    |  @   a   GET_YMID_R_W%I    �  @   a   GET_YMID_R_W%J    �  @   a   GET_YMID_R_W%L    <   e       GET_YSIN    �   @   a   GET_YSIN%I    �   @   a   GET_YSIN%J    !!  @   a   GET_YSIN%L    a!  u       GET_XOFFSET $   �!  @      GET_XOFFSET%PRESENT #   "  @   a   GET_XOFFSET%GLOBAL    V"  u       GET_YOFFSET $   �"  @      GET_YOFFSET%PRESENT #   #  @   a   GET_YOFFSET%GLOBAL    K#  �       INIT_GRID $   �#  @   a   INIT_GRID%AM_I_ROOT    $  @   a   INIT_GRID%IM    ^$  @   a   INIT_GRID%JM    �$  @   a   INIT_GRID%LM    �$  @   a   INIT_GRID%RC "   %  P       ITS_A_NESTED_GRID    n%  V       SET_XOFFSET %   �%  @   a   SET_XOFFSET%X_OFFSET    &  V       SET_YOFFSET %   Z&  @   a   SET_YOFFSET%Y_OFFSET 