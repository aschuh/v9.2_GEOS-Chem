	  �Q  �   k820309    V          13.1        lOaT                                                                                                           
       time_mod.F TIME_MOD       l       SET_CURRENT_TIME SET_BEGIN_TIME SET_END_TIME SET_NDIAGTIME SET_DIAGB SET_DIAGE SET_TIMESTEPS SET_CT_CHEM SET_CT_CONV SET_CT_DYN SET_CT_EMIS SET_CT_DIAG SET_CT_A1 SET_CT_A3 SET_CT_A6 SET_CT_I3 SET_CT_I6 SET_CT_XTRA SET_ELAPSED_MIN GET_JD GET_ELAPSED_MIN GET_ELAPSED_SEC GET_NYMDB GET_NHMSB GET_NYMDE GET_NHMSE GET_NYMD GET_NHMS GET_NDIAGTIME GET_TIME_AHEAD GET_MONTH GET_DAY GET_YEAR GET_HOUR GET_MINUTE GET_SECOND GET_DAY_OF_YEAR GET_DAY_OF_WEEK GET_DAY_OF_WEEK_LT GET_GMT GET_TAU GET_TAUB GET_TAUE GET_DIAGB GET_DIAGE GET_LOCALTIME GET_SEASON GET_TS_CHEM GET_TS_CONV GET_TS_DIAG GET_TS_DYN GET_TS_EMIS GET_TS_UNIT GET_CT_CHEM GET_CT_CONV GET_CT_DYN GET_CT_EMIS GET_CT_A1 GET_CT_A3 GET_CT_A6 GET_CT_I3 GET_CT_I6 GET_CT_XTRA GET_CT_DIAG GET_A1_TIME GET_A3_TIME GET_A6_TIME GET_I3_TIME GET_I6_TIME GET_FIRST_A1_TIME GET_FIRST_A3_TIME GET_FIRST_A6_TIME GET_FIRST_I3_TIME GET_FIRST_I6_TIME ITS_TIME_FOR_CHEM ITS_TIME_FOR_CONV ITS_TIME_FOR_DYN ITS_TIME_FOR_EMIS ITS_TIME_FOR_UNIT ITS_TIME_FOR_DIAG ITS_TIME_FOR_A1 ITS_TIME_FOR_A3 ITS_TIME_FOR_A6 ITS_TIME_FOR_I3 ITS_TIME_FOR_I6 ITS_TIME_FOR_UNZIP ITS_TIME_FOR_DEL ITS_TIME_FOR_EXIT ITS_TIME_FOR_BPCH ITS_A_LEAPYEAR ITS_A_NEW_YEAR ITS_A_NEW_MONTH ITS_MIDMONTH ITS_A_NEW_DAY ITS_A_NEW_HOUR ITS_A_NEW_SEASON PRINT_CURRENT_TIME TIMESTAMP_STRING YMD_EXTRACT EXPAND_DATE SYSTEM_DATE_TIME SYSTEM_TIMESTAMP TIMESTAMP_DIAG GET_NYMD_DIAG GET_HG2_DIAG SET_HG2_DIAG SET_HISTYR GET_HISTYR #         @                                                       #SET_CURRENT_TIME%SIGN    #SET_CURRENT_TIME%NINT    #SET_CURRENT_TIME%INT    #SET_CURRENT_TIME%MOD    #SET_CURRENT_TIME%DBLE                  @                                 SIGN               @                                 NINT               @                                 INT               @                                 MOD               @                                 DBLE #         @                                                      #SET_BEGIN_TIME%DBLE    #THISNYMDB 	   #THISNHMSB 
                 @                                 DBLE           
   @                              	                     
   @                              
           #         @                                                      #SET_END_TIME%DBLE    #THISNYMDE    #THISNHMSE                  @                                 DBLE           
   @                                                   
   @                                         #         @                                                       #THIS_NDIAGTIME              
   @                                         #         @                                                       #THISDIAGB              
   @                                  
      #         @                                                       #THISDIAGE              
   @                                  
      #         @                                                       #AM_I_ROOT    #CHEMISTRY    #CONVECTION    #DYNAMICS    #EMISSION    #UNIT_CONV    #DIAGNOS              
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                         #         @                                                      #SET_CT_CHEM%PRESENT    #INCREMENT    #RESET                   @                                 PRESENT           
 @@                                                   
 @@                                          #         @                                   !                   #SET_CT_CONV%PRESENT "   #INCREMENT #   #RESET $                 @                            "     PRESENT           
 @@                              #                     
 @@                              $           #         @                                   %                   #SET_CT_DYN%PRESENT &   #INCREMENT '   #RESET (                 @                            &     PRESENT           
 @@                              '                     
 @@                              (           #         @                                   )                   #SET_CT_EMIS%PRESENT *   #INCREMENT +   #RESET ,                 @                            *     PRESENT           
 @@                              +                     
 @@                              ,           #         @                                   -                   #SET_CT_DIAG%PRESENT .   #INCREMENT /   #RESET 0                 @                            .     PRESENT           
 @@                              /                     
 @@                              0           #         @                                   1                   #SET_CT_A1%PRESENT 2   #INCREMENT 3   #RESET 4                 @                            2     PRESENT           
 @@                              3                     
 @@                              4           #         @                                   5                   #SET_CT_A3%PRESENT 6   #INCREMENT 7   #RESET 8                 @                            6     PRESENT           
 @@                              7                     
 @@                              8           #         @                                   9                   #SET_CT_A6%PRESENT :   #INCREMENT ;   #RESET <                 @                            :     PRESENT           
 @@                              ;                     
 @@                              <           #         @                                   =                   #SET_CT_I3%PRESENT >   #INCREMENT ?   #RESET @                 @                            >     PRESENT           
 @@                              ?                     
 @@                              @           #         @                                   A                   #SET_CT_I6%PRESENT B   #INCREMENT C   #RESET D                 @                            B     PRESENT           
 @@                              C                     
 @@                              D           #         @                                   E                   #SET_CT_XTRA%PRESENT F   #INCREMENT G   #RESET H                 @                            F     PRESENT           
 @@                              G                     
 @@                              H           #         @                                   I                     %         @                                J                   
       #GET_JD%DBLE K   #THISNYMD L   #THISNHMS M                 @                            K     DBLE           
  @@                              L                     
  @@                              M           %         @                                 N                            %         @                                 O                            %         @                                 P                            %         @                                 Q                            %         @                                 R                            %         @                                 S                            %         @                                T                            %         @                                 U                            %         @                                 V                            (         `                                W                                       #N_MINS X   p          p            p                                    
   @                              X           %         @                                 Y                            %         @                                 Z                            %         @                                 [                            %         @                                 \                            %         @                                 ]                            %         @                                 ^                            %         @                                 _                            %         @                                 `                            %         @                                 a                          #GET_DAY_OF_WEEK_LT%INT b   #GET_DAY_OF_WEEK_LT%MOD c   #GET_DAY_OF_WEEK_LT%DBLE d   #I e   #J f   #L g                 @                            b     INT               @                            c     MOD               @                            d     DBLE           
  @@                              e                     
  @@                              f                     
  @@                              g           %         @                                h                     
       %         @                                 i                     
       %         @                                 j                     
       %         @                                 k                     
       %         @                                 l                            %         @                                 m                            %         @                                 n                   
       #GET_LOCALTIME%PRESENT o   #I p   #J q   #L r   #GMT s                 @                            o     PRESENT           
  @@                              p                     
  @@                              q                     
  @@                              r                     
 @@                             s     
      %         @                                 t                            %         @                                 u                            %         @                                 v                            %         @                                 w                            %         @                                 x                            %         @                                 y                            %         @                                 z                            %         @                                 {                            %         @                                 |                            %         @                                 }                            %         @                                 ~                            %         @                                                             %         @                                 �                            %         @                                 �                            %         @                                 �                            %         @                                 �                            %         @                                 �                            %         @                                 �                            (         `                                �                                        p          p            p                          (         `                                 �                                        p          p            p                          (         `                                 �                                        p          p            p                          (         `                                 �                                       #GET_I3_TIME%MOD �   p          p            p                                        @                            �     MOD (         `                                 �                                       #GET_I6_TIME%MOD �   p          p            p                                        @                            �     MOD (         `                                 �                                        p          p            p                          (         `                                 �                                       #GET_FIRST_A3_TIME%MOD �   p          p            p                                        @                            �     MOD (         `                                 �                   	                    #GET_FIRST_A6_TIME%MOD �   p          p            p                                        @                            �     MOD (         `                                 �                   
                    #GET_FIRST_I3_TIME%MOD �   p          p            p                                        @                            �     MOD (         `                                 �                                       #GET_FIRST_I6_TIME%MOD �   p          p            p                                        @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_CHEM%MAX �   #ITS_TIME_FOR_CHEM%MOD �                 @                            �     MAX               @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_CONV%MOD �                 @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_DYN%MOD �                 @                            �     MOD %         @                                �                           #ITS_TIME_FOR_EMIS%MAX �   #ITS_TIME_FOR_EMIS%MOD �                 @                            �     MAX               @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_UNIT%MOD �                 @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_DIAG%MOD �                 @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_A1%MOD �                 @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_A3%MOD �                 @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_A6%MOD �                 @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_I3%MOD �                 @                            �     MOD %         @                                 �                           #ITS_TIME_FOR_I6%MOD �                 @                            �     MOD %         @                                 �                            %         @                                 �                            %         @                                 �                            %         @                                 �                            %         @                                �                          #ITS_A_LEAPYEAR%PRESENT �   #ITS_A_LEAPYEAR%MOD �   #YEAR_IN �   #FORCE �                 @                            �     PRESENT               @                            �     MOD           
 @@                              �                     
 @@                              �           %         @                                 �                          #ITS_A_NEW_YEAR%PRESENT �   #ITS_A_NEW_YEAR%DBLE �   #NO_CCTS �                 @                            �     PRESENT               @                            �     DBLE            @@                              �            %         @                                 �                          #ITS_A_NEW_MONTH%PRESENT �   #ITS_A_NEW_MONTH%DBLE �   #NO_CCTS �                 @                            �     PRESENT               @                            �     DBLE            @@                              �            %         @                                 �                            %         @                                �                          #ITS_A_NEW_DAY%PRESENT �   #ITS_A_NEW_DAY%DBLE �   #NO_CCTS �                 @                            �     PRESENT               @                            �     DBLE            @@                              �            %         @                                 �                           #ITS_A_NEW_HOUR%MOD �                 @                            �     MOD %         @                                 �                            #         @                                   �                    #PRINT_CURRENT_TIME%REAL �                 @                            �     REAL $         @                                �                          #TIMESTAMP_STRING%PRESENT �   #YYYYMMDD �   #HHMMSS �                         @                            �     PRESENT           
 @@                              �                     
 @@                              �           #         @                                  �                   #YMD_EXTRACT%INT �   #YMD_EXTRACT%DBLE �   #NYMD �   #Y �   #M �   #D �                 @                            �     INT               @                            �     DBLE           
  @@                              �                     D @@                              �                      D @@                              �                      D  @                              �            #         @                                   �                    #FILENAME �   #YYYYMMDD �   #HHMMSS �             
D @@                             �                     1           
  @@                              �                     
  @@                              �           #         @                                  �                    #SYS_NYMD �   #SYS_NHMS �                                                 D  @                              �                      D  @                              �            $         @                                 �                                    #         @                                   �                     %         @                                 �                            %         @                                 �                            #         @                                   �                   #SET_HG2_DIAG%PRESENT �   #INCREMENT �   #RESET �                 @                            �     PRESENT           
 @@                              �                     
 @@                              �           #         @                                   �                    #YEARIN �             
   @                              �           %         @                                 �                               �         fn#fn    �   �  b   uapp(TIME_MOD !   Y  �       SET_CURRENT_TIME &   &  =      SET_CURRENT_TIME%SIGN &   c  =      SET_CURRENT_TIME%NINT %   �  <      SET_CURRENT_TIME%INT %   �  <      SET_CURRENT_TIME%MOD &     =      SET_CURRENT_TIME%DBLE    U         SET_BEGIN_TIME $   �  =      SET_BEGIN_TIME%DBLE )   	  @   a   SET_BEGIN_TIME%THISNYMDB )   Q	  @   a   SET_BEGIN_TIME%THISNHMSB    �	  }       SET_END_TIME "   
  =      SET_END_TIME%DBLE '   K
  @   a   SET_END_TIME%THISNYMDE '   �
  @   a   SET_END_TIME%THISNHMSE    �
  \       SET_NDIAGTIME -   '  @   a   SET_NDIAGTIME%THIS_NDIAGTIME    g  W       SET_DIAGB $   �  @   a   SET_DIAGB%THISDIAGB    �  W       SET_DIAGE $   U  @   a   SET_DIAGE%THISDIAGE    �  �       SET_TIMESTEPS (   C  @   a   SET_TIMESTEPS%AM_I_ROOT (   �  @   a   SET_TIMESTEPS%CHEMISTRY )   �  @   a   SET_TIMESTEPS%CONVECTION '     @   a   SET_TIMESTEPS%DYNAMICS '   C  @   a   SET_TIMESTEPS%EMISSION (   �  @   a   SET_TIMESTEPS%UNIT_CONV &   �  @   a   SET_TIMESTEPS%DIAGNOS      {       SET_CT_CHEM $   ~  @      SET_CT_CHEM%PRESENT &   �  @   a   SET_CT_CHEM%INCREMENT "   �  @   a   SET_CT_CHEM%RESET    >  {       SET_CT_CONV $   �  @      SET_CT_CONV%PRESENT &   �  @   a   SET_CT_CONV%INCREMENT "   9  @   a   SET_CT_CONV%RESET    y  z       SET_CT_DYN #   �  @      SET_CT_DYN%PRESENT %   3  @   a   SET_CT_DYN%INCREMENT !   s  @   a   SET_CT_DYN%RESET    �  {       SET_CT_EMIS $   .  @      SET_CT_EMIS%PRESENT &   n  @   a   SET_CT_EMIS%INCREMENT "   �  @   a   SET_CT_EMIS%RESET    �  {       SET_CT_DIAG $   i  @      SET_CT_DIAG%PRESENT &   �  @   a   SET_CT_DIAG%INCREMENT "   �  @   a   SET_CT_DIAG%RESET    )  y       SET_CT_A1 "   �  @      SET_CT_A1%PRESENT $   �  @   a   SET_CT_A1%INCREMENT     "  @   a   SET_CT_A1%RESET    b  y       SET_CT_A3 "   �  @      SET_CT_A3%PRESENT $     @   a   SET_CT_A3%INCREMENT     [  @   a   SET_CT_A3%RESET    �  y       SET_CT_A6 "     @      SET_CT_A6%PRESENT $   T  @   a   SET_CT_A6%INCREMENT     �  @   a   SET_CT_A6%RESET    �  y       SET_CT_I3 "   M  @      SET_CT_I3%PRESENT $   �  @   a   SET_CT_I3%INCREMENT     �  @   a   SET_CT_I3%RESET      y       SET_CT_I6 "   �  @      SET_CT_I6%PRESENT $   �  @   a   SET_CT_I6%INCREMENT       @   a   SET_CT_I6%RESET    F  {       SET_CT_XTRA $   �  @      SET_CT_XTRA%PRESENT &     @   a   SET_CT_XTRA%INCREMENT "   A  @   a   SET_CT_XTRA%RESET     �  H       SET_ELAPSED_MIN    �  }       GET_JD    F  =      GET_JD%DBLE     �  @   a   GET_JD%THISNYMD     �  @   a   GET_JD%THISNHMS       P       GET_ELAPSED_MIN     S  P       GET_ELAPSED_SEC    �  P       GET_NYMDB    �  P       GET_NHMSB    C  P       GET_NYMDE    �  P       GET_NHMSE    �  P       GET_NYMD    3   P       GET_NHMS    �   P       GET_NDIAGTIME    �   �       GET_TIME_AHEAD &   �!  @   a   GET_TIME_AHEAD%N_MINS    �!  P       GET_MONTH    "  P       GET_DAY    c"  P       GET_YEAR    �"  P       GET_HOUR    #  P       GET_MINUTE    S#  P       GET_SECOND     �#  P       GET_DAY_OF_YEAR     �#  P       GET_DAY_OF_WEEK #   C$  �       GET_DAY_OF_WEEK_LT '   �$  <      GET_DAY_OF_WEEK_LT%INT '   9%  <      GET_DAY_OF_WEEK_LT%MOD (   u%  =      GET_DAY_OF_WEEK_LT%DBLE %   �%  @   a   GET_DAY_OF_WEEK_LT%I %   �%  @   a   GET_DAY_OF_WEEK_LT%J %   2&  @   a   GET_DAY_OF_WEEK_LT%L    r&  P       GET_GMT    �&  P       GET_TAU    '  P       GET_TAUB    b'  P       GET_TAUE    �'  P       GET_DIAGB    (  P       GET_DIAGE    R(  �       GET_LOCALTIME &   �(  @      GET_LOCALTIME%PRESENT     )  @   a   GET_LOCALTIME%I     [)  @   a   GET_LOCALTIME%J     �)  @   a   GET_LOCALTIME%L "   �)  @   a   GET_LOCALTIME%GMT    *  P       GET_SEASON    k*  P       GET_TS_CHEM    �*  P       GET_TS_CONV    +  P       GET_TS_DIAG    [+  P       GET_TS_DYN    �+  P       GET_TS_EMIS    �+  P       GET_TS_UNIT    K,  P       GET_CT_CHEM    �,  P       GET_CT_CONV    �,  P       GET_CT_DYN    ;-  P       GET_CT_EMIS    �-  P       GET_CT_A1    �-  P       GET_CT_A3    +.  P       GET_CT_A6    {.  P       GET_CT_I3    �.  P       GET_CT_I6    /  P       GET_CT_XTRA    k/  P       GET_CT_DIAG    �/  �       GET_A1_TIME    _0  �       GET_A3_TIME    1  �       GET_A6_TIME    �1  �       GET_I3_TIME     `2  <      GET_I3_TIME%MOD    �2  �       GET_I6_TIME     U3  <      GET_I6_TIME%MOD "   �3  �       GET_FIRST_A1_TIME "   54  �       GET_FIRST_A3_TIME &   �4  <      GET_FIRST_A3_TIME%MOD "   05  �       GET_FIRST_A6_TIME &   �5  <      GET_FIRST_A6_TIME%MOD "   +6  �       GET_FIRST_I3_TIME &   �6  <      GET_FIRST_I3_TIME%MOD "   &7  �       GET_FIRST_I6_TIME &   �7  <      GET_FIRST_I6_TIME%MOD "   !8  �       ITS_TIME_FOR_CHEM &   �8  <      ITS_TIME_FOR_CHEM%MAX &   �8  <      ITS_TIME_FOR_CHEM%MOD "   9  k       ITS_TIME_FOR_CONV &   �9  <      ITS_TIME_FOR_CONV%MOD !   �9  j       ITS_TIME_FOR_DYN %   0:  <      ITS_TIME_FOR_DYN%MOD "   l:  �       ITS_TIME_FOR_EMIS &   �:  <      ITS_TIME_FOR_EMIS%MAX &   .;  <      ITS_TIME_FOR_EMIS%MOD "   j;  k       ITS_TIME_FOR_UNIT &   �;  <      ITS_TIME_FOR_UNIT%MOD "   <  k       ITS_TIME_FOR_DIAG &   |<  <      ITS_TIME_FOR_DIAG%MOD     �<  i       ITS_TIME_FOR_A1 $   !=  <      ITS_TIME_FOR_A1%MOD     ]=  i       ITS_TIME_FOR_A3 $   �=  <      ITS_TIME_FOR_A3%MOD     >  i       ITS_TIME_FOR_A6 $   k>  <      ITS_TIME_FOR_A6%MOD     �>  i       ITS_TIME_FOR_I3 $   ?  <      ITS_TIME_FOR_I3%MOD     L?  i       ITS_TIME_FOR_I6 $   �?  <      ITS_TIME_FOR_I6%MOD #   �?  P       ITS_TIME_FOR_UNZIP !   A@  P       ITS_TIME_FOR_DEL "   �@  P       ITS_TIME_FOR_EXIT "   �@  P       ITS_TIME_FOR_BPCH    1A  �       ITS_A_LEAPYEAR '   �A  @      ITS_A_LEAPYEAR%PRESENT #   B  <      ITS_A_LEAPYEAR%MOD '   IB  @   a   ITS_A_LEAPYEAR%YEAR_IN %   �B  @   a   ITS_A_LEAPYEAR%FORCE    �B  �       ITS_A_NEW_YEAR '   [C  @      ITS_A_NEW_YEAR%PRESENT $   �C  =      ITS_A_NEW_YEAR%DBLE '   �C  @   a   ITS_A_NEW_YEAR%NO_CCTS     D  �       ITS_A_NEW_MONTH (   �D  @      ITS_A_NEW_MONTH%PRESENT %   �D  =      ITS_A_NEW_MONTH%DBLE (   )E  @   a   ITS_A_NEW_MONTH%NO_CCTS    iE  P       ITS_MIDMONTH    �E  �       ITS_A_NEW_DAY &   IF  @      ITS_A_NEW_DAY%PRESENT #   �F  =      ITS_A_NEW_DAY%DBLE &   �F  @   a   ITS_A_NEW_DAY%NO_CCTS    G  h       ITS_A_NEW_HOUR #   nG  <      ITS_A_NEW_HOUR%MOD !   �G  P       ITS_A_NEW_SEASON #   �G  e       PRINT_CURRENT_TIME (   _H  =      PRINT_CURRENT_TIME%REAL !   �H  �       TIMESTAMP_STRING )   ,I  @      TIMESTAMP_STRING%PRESENT *   lI  @   a   TIMESTAMP_STRING%YYYYMMDD (   �I  @   a   TIMESTAMP_STRING%HHMMSS    �I  �       YMD_EXTRACT     ~J  <      YMD_EXTRACT%INT !   �J  =      YMD_EXTRACT%DBLE !   �J  @   a   YMD_EXTRACT%NYMD    7K  @   a   YMD_EXTRACT%Y    wK  @   a   YMD_EXTRACT%M    �K  @   a   YMD_EXTRACT%D    �K  p       EXPAND_DATE %   gL  L   a   EXPAND_DATE%FILENAME %   �L  @   a   EXPAND_DATE%YYYYMMDD #   �L  @   a   EXPAND_DATE%HHMMSS !   3M  �       SYSTEM_DATE_TIME *   �M  @   a   SYSTEM_DATE_TIME%SYS_NYMD *   �M  @   a   SYSTEM_DATE_TIME%SYS_NHMS !   ;N  X       SYSTEM_TIMESTAMP    �N  H       TIMESTAMP_DIAG    �N  P       GET_NYMD_DIAG    +O  P       GET_HG2_DIAG    {O  |       SET_HG2_DIAG %   �O  @      SET_HG2_DIAG%PRESENT '   7P  @   a   SET_HG2_DIAG%INCREMENT #   wP  @   a   SET_HG2_DIAG%RESET    �P  T       SET_HISTYR "   Q  @   a   SET_HISTYR%YEARIN    KQ  P       GET_HISTYR 