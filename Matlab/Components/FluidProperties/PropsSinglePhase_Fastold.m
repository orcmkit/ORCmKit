function [mu, Pr, k, mu_rat] = PropsSinglePhase_Fast(T_mean, P_mean, T_wall, T_sat_mean, zone, fluid)


if strcmp(fluid(1:3), 'ICP')
    k =  PropsSI_ICP('L', 'T', T_mean, 'P', P_mean, fluid);
    mu = PropsSI_ICP('V', 'T', T_mean, 'P', P_mean, fluid);
    cp = PropsSI_ICP('C', 'T', T_mean, 'P', P_mean, fluid);
    Pr = cp*mu/k;
    mu_wall = PropsSI_ICP('V', 'T', T_wall, 'P', P_mean, fluid);
    mu_rat = mu/mu_wall;
else
    if (strcmp(fluid, 'R245fa'))
        if (strcmp(zone, 'liq')) && (T_mean >= T_sat_mean-5e-2) %--> DONE
            T_mean = max(253.15, min(T_mean, 423.15));                        
            mean = 0; std = 1; mu = (-1.21860675033688e-14)*((T_mean-mean)/std)^5+(2.17892205367452e-11)*((T_mean-mean)/std)^4+(-1.56320390614916e-08)*((T_mean-mean)/std)^3+(5.63657180203827e-06)*((T_mean-mean)/std)^2+(-0.00102596558424522)*((T_mean-mean)/std)+(0.0760785978040368);            
            mean = 3.801453576303341e+02; std = 38.247859808756843; k  = (0)*exp(-((((T_mean-mean)/std)-(-7.59170750439347))/(0.651117339142076))^2)+(0.109163296670942)*exp(-((((T_mean-mean)/std)-(-4.20344872532284))/(2.34246199427738))^2)+(-0.0769406610602091)*exp(-((((T_mean-mean)/std)-(-1.10772542059559))/(1.30809883808847))^2)+(0.1321402268234)*exp(-((((T_mean-mean)/std)-(-1.08598655946568))/(1.50327347131777))^2)+(0.0405221286285368)*exp(-((((T_mean-mean)/std)-(1.33577807704584))/(1.58688382667343))^2)+(1293783626.20023)*exp(-((((T_mean-mean)/std)-(4.05348169135316))/(0.563304682670529))^2);            
            mean = 3.748225461986881e+02; std = 37.219304126609593; cp = (5.08722856861344e+15)*exp(-((((T_mean-mean)/std)-(7.70351575334858))/(1.19675531435243))^2)+(166.867912055228)*exp(-((((T_mean-mean)/std)-(1.20928277711077))/(0.409573169659633))^2)+(490519.551970258)*exp(-((((T_mean-mean)/std)-(10.9770110099433))/(3.79811865409817))^2)+(-2.36802627255662)*exp(-((((T_mean-mean)/std)-(0.0794360924739605))/(0.0852293679390625))^2)+(0)*exp(-((((T_mean-mean)/std)-(-10.954160642599))/(0.483931776153539))^2)+(5759.80654330883)*exp(-((((T_mean-mean)/std)-(45.2285480580039))/(38.9147051615508))^2);
            Pr = cp*mu/k;
            mu_rat = 1;
            
        elseif (strcmp(zone, 'liq')) && (T_mean < T_sat_mean-5e-2)   %--> DONE
            T_mean = max(253.15, min(T_mean, 423.15)); P_mean = max(0.1e5, min(P_mean, 36e5));  T_wall = max(253.15, min(T_wall, 423.15));                       
            mean_T = 3.162660832157676e+02; std_T = 43.137399521369431;  mean_P = 1.804999999999939e+06; std_P = 1.046811851712497e+06; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; mu = (0.000333171199094998)+(-0.000171569020839364)*(T_mean_bis)+(5.86291135254676e-06)*(P_mean_bis)+(5.1929796694257e-05)*(T_mean_bis)^2+(-8.68696115745206e-07)*(T_mean_bis)*(P_mean_bis)+(2.36191751292916e-08)*(P_mean_bis)^2+(-1.97116873362024e-05)*(T_mean_bis)^3+(5.59305116003848e-07)*(T_mean_bis)^2*(P_mean_bis)+(-1.16440486750793e-07)*(T_mean_bis)*(P_mean_bis)^2+(7.51596530135219e-06)*(T_mean_bis)^4+(-4.58208899040354e-07)*(T_mean_bis)^3*(P_mean_bis)+(-1.28579794615051e-07)*(T_mean_bis)^2*(P_mean_bis)^2+(-1.83974283613259e-06)*(T_mean_bis)^5+(5.05336639965252e-07)*(T_mean_bis)^4*(P_mean_bis)+(5.0578520156275e-08)*(T_mean_bis)^3*(P_mean_bis)^2;
            mean_T = 3.164871827726678e+02; std_T = 42.994186075806361;  mean_P = 1.809999999999994e+06; std_P = 1.043895941262018e+06; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; k = (0.0854472601116622)+(-0.0132841940252005)*(T_mean_bis)+(0.00061296884420845)*(P_mean_bis)+(-9.99484363934673e-05)*(T_mean_bis)^2+(0.000208642029420873)*(T_mean_bis)*(P_mean_bis)+(-3.33995693551934e-06)*(P_mean_bis)^2;
            if P_mean <=20e5
                mean_T = 3.039858691039196e+02; std_T = 35.022443669207711;  mean_P = 1010000; std_P = 5.773502691896258e+05; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; cp = (1327.94914261602)+(88.3075001950569)*(T_mean_bis)+(-2.06135367407864)*(P_mean_bis)+(11.581039403172)*(T_mean_bis)^2+(-1.28982763559531)*(T_mean_bis)*(P_mean_bis)+(-0.0620540373536512)*(P_mean_bis)^2+(-1.11148954053076)*(T_mean_bis)^3+(-0.667977436837543)*(T_mean_bis)^2*(P_mean_bis)+(-0.165069555263178)*(T_mean_bis)*(P_mean_bis)^2+(0.00984012053428237)*(P_mean_bis)^3+(1.04237267752328)*(T_mean_bis)^4+(-0.789628025087167)*(T_mean_bis)^3*(P_mean_bis)+(0.0600253111616615)*(T_mean_bis)^2*(P_mean_bis)^2+(0.034945872637925)*(T_mean_bis)*(P_mean_bis)^3+(0.0294440240379889)*(P_mean_bis)^4+(1.64418763436779)*(T_mean_bis)^5+(-0.416565894140683)*(T_mean_bis)^4*(P_mean_bis)+(0.125834782536733)*(T_mean_bis)^3*(P_mean_bis)^2+(0.0623503865887964)*(T_mean_bis)^2*(P_mean_bis)^3+(0.0100636311832893)*(T_mean_bis)*(P_mean_bis)^4+(-0.0209798751778806)*(P_mean_bis)^5;
            else
                mean_T = 3.322338429969005e+02; std_T =  46.419407482435048;  mean_P =  2.800000000000033e+06; std_P =  4.665456720724282e+05; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; cp = (1397.98884805449)+(137.472521600461)*(T_mean_bis)+(-1.91516409643155)*(P_mean_bis)+(4.75711476000692)*(T_mean_bis)^2+(-4.27010318564583)*(T_mean_bis)*(P_mean_bis)+(-0.0356048138065627)*(P_mean_bis)^2+(-2.06991946101793)*(T_mean_bis)^3+(-7.75019850622227)*(T_mean_bis)^2*(P_mean_bis)+(-0.146228283476242)*(T_mean_bis)*(P_mean_bis)^2+(-0.0539687200643268)*(P_mean_bis)^3+(28.5218958742178)*(T_mean_bis)^4+(-3.69621327638122)*(T_mean_bis)^3*(P_mean_bis)+(0.0842944280572547)*(T_mean_bis)^2*(P_mean_bis)^2+(0.0939875013662841)*(T_mean_bis)*(P_mean_bis)^3+(0.0371094376707812)*(P_mean_bis)^4+(15.2524843062994)*(T_mean_bis)^5+(-0.00414426226109706)*(T_mean_bis)^4*(P_mean_bis)+(0.108014654870543)*(T_mean_bis)^3*(P_mean_bis)^2+(0.107274821766002)*(T_mean_bis)^2*(P_mean_bis)^3+(0.0498772212270696)*(T_mean_bis)*(P_mean_bis)^4+(0.00719439141040417)*(P_mean_bis)^5;
            end
            Pr = cp*mu/k;
            mean_T = 3.162660832157676e+02; std_T = 43.137399521369431;  mean_P = 1.804999999999939e+06; std_P = 1.046811851712497e+06; T_mean_bis = (T_wall - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; mu_wall = (0.000333171199094998)+(-0.000171569020839364)*(T_mean_bis)+(5.86291135254676e-06)*(P_mean_bis)+(5.1929796694257e-05)*(T_mean_bis)^2+(-8.68696115745206e-07)*(T_mean_bis)*(P_mean_bis)+(2.36191751292916e-08)*(P_mean_bis)^2+(-1.97116873362024e-05)*(T_mean_bis)^3+(5.59305116003848e-07)*(T_mean_bis)^2*(P_mean_bis)+(-1.16440486750793e-07)*(T_mean_bis)*(P_mean_bis)^2+(7.51596530135219e-06)*(T_mean_bis)^4+(-4.58208899040354e-07)*(T_mean_bis)^3*(P_mean_bis)+(-1.28579794615051e-07)*(T_mean_bis)^2*(P_mean_bis)^2+(-1.83974283613259e-06)*(T_mean_bis)^5+(5.05336639965252e-07)*(T_mean_bis)^4*(P_mean_bis)+(5.0578520156275e-08)*(T_mean_bis)^3*(P_mean_bis)^2;   
            mu_rat = mu/mu_wall;

        elseif (strcmp(zone, 'vap') || strcmp(zone, 'vap_wet') || strcmp(zone, 'tp_dryout') ) && (T_mean <= T_sat_mean +5e-2) %--> DONE
            T_mean = max(253.15, min(T_mean, 423.15));                        
            mean = 3.801453576303341e+02; std = 38.247859808756843; mu = (3.80375516197253e-08)*((T_mean-mean)/std)^8+(3.82704211085208e-07)*((T_mean-mean)/std)^7+(1.33104510066901e-06)*((T_mean-mean)/std)^6+(1.64813604083368e-06)*((T_mean-mean)/std)^5+(-1.67207038516526e-07)*((T_mean-mean)/std)^4+(-9.14356388317895e-07)*((T_mean-mean)/std)^3+(1.21659693562957e-06)*((T_mean-mean)/std)^2+(2.89121162408815e-06)*((T_mean-mean)/std)+(1.33545184937567e-05);
            mean = 3.789778721533120e+02; std = 38.023629238075770; cp =  (1.20949478768067e+16)*exp(-((((T_mean-mean)/std)-(3.01996185889906))/(0.337174427496613))^2)+(2182207898902.77)*exp(-((((T_mean-mean)/std)-(6.4681966584702))/(1.17638815482174))^2)+(13.7614693833771)*exp(-((((T_mean-mean)/std)-(0.771111960694984))/(0.097646393467277))^2)+(-7.7177041017477)*exp(-((((T_mean-mean)/std)-(0.58956921717628))/(0.0325884030990685))^2)+(1702.62546483762)*exp(-((((T_mean-mean)/std)-(1.84750833796286))/(1.05752033508848))^2)+(6288410.86494768)*exp(-((((T_mean-mean)/std)-(113.035430082245))/(38.7138677040909))^2);          
            mean = 3.801453576303341e+02; std = 38.247859808756843; k  = (8809.78480089158)*exp(-((((T_mean-mean)/std)-(1.85180384630831))/(0.175020541150224))^2)+(0.149343590838944)*exp(-((((T_mean-mean)/std)-(1.90490653426044))/(0.48593795772417))^2)+(0.0062664265487332)*exp(-((((T_mean-mean)/std)-(1.1225364883482))/(0.529461525267281))^2)+(0.00993647246617092)*exp(-((((T_mean-mean)/std)-(1.42839159349335))/(1.49963297811562))^2)+(0)*exp(-((((T_mean-mean)/std)-(0.199552371607423))/(0.000459173301811569))^2)+(0.0189415487287949)*exp(-((((T_mean-mean)/std)-(0.766680937523611))/(5.25502326542243))^2);
            Pr = cp*mu/k;
            mu_rat = 1;
            
        elseif (strcmp(zone, 'vap') || strcmp(zone, 'vap_wet') || strcmp(zone, 'tp_dryout') ) && (T_mean > T_sat_mean +5e-2) %--> DONE
            T_mean = max(253.15, min(T_mean, 423.15)); P_mean = max(0.1e5, min(P_mean, 36e5)); T_wall = max(253.15, min(T_wall, 423.15));                       
            mean_T = 4.046040468982957e+02; std_T = 24.924959340389851;  mean_P = 1.848966986417469e+06; std_P = 1.023850577380752e+06; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; mu = (1.44569298567578e-05)+(8.04897653766163e-07)*(T_mean_bis)+(1.2046133534143e-06)*(P_mean_bis)+(3.28204169912859e-07)*(T_mean_bis)^2+(-8.13776321013716e-07)*(T_mean_bis)*(P_mean_bis)+(1.02088562667795e-06)*(P_mean_bis)^2+(-1.05166329365479e-07)*(T_mean_bis)^3+(1.25394383389519e-06)*(T_mean_bis)^2*(P_mean_bis)+(-2.48827554255973e-06)*(T_mean_bis)*(P_mean_bis)^2+(1.19193622360903e-06)*(P_mean_bis)^3+(-9.33125361157745e-09)*(T_mean_bis)^4+(-1.48694274071809e-07)*(T_mean_bis)^3*(P_mean_bis)+(1.43611038860212e-06)*(T_mean_bis)^2*(P_mean_bis)^2+(-2.76907108924106e-06)*(T_mean_bis)*(P_mean_bis)^3+(1.34110616478549e-06)*(P_mean_bis)^4+(7.14653026436247e-10)*(T_mean_bis)^5+(-9.89795305435064e-09)*(T_mean_bis)^4*(P_mean_bis)+(-4.37456744704995e-08)*(T_mean_bis)^3*(P_mean_bis)^2+(4.8157748815481e-07)*(T_mean_bis)^2*(P_mean_bis)^3+(-9.42084427674813e-07)*(T_mean_bis)*(P_mean_bis)^4+(4.76457033599083e-07)*(P_mean_bis)^5;        
            mean_T = 4.046040468982957e+02; std_T = 24.924959340389851;  mean_P = 1.848966986417469e+06; std_P = 1.023850577380752e+06; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; k = (0.0247682795644727)+(0.000649716344521826)*(T_mean_bis)+(0.00403608882569339)*(P_mean_bis)+(0.00123433463255559)*(T_mean_bis)^2+(-0.00459931517278111)*(T_mean_bis)*(P_mean_bis)+(0.00389901634029654)*(P_mean_bis)^2+(-0.00060565169082321)*(T_mean_bis)^3+(0.0046009699624889)*(T_mean_bis)^2*(P_mean_bis)+(-0.00846909177658402)*(T_mean_bis)*(P_mean_bis)^2+(0.0040160122311898)*(P_mean_bis)^3+(6.0920410024453e-05)*(T_mean_bis)^4+(-0.00116694120396312)*(T_mean_bis)^3*(P_mean_bis)+(0.00535310527893279)*(T_mean_bis)^2*(P_mean_bis)^2+(-0.00763836924448965)*(T_mean_bis)*(P_mean_bis)^3+(0.00305067420242232)*(P_mean_bis)^4+(-4.03219418189149e-07)*(T_mean_bis)^5+(4.24909343714003e-05)*(T_mean_bis)^4*(P_mean_bis)+(-0.000511813924430031)*(T_mean_bis)^3*(P_mean_bis)^2+(0.00188041165149096)*(T_mean_bis)^2*(P_mean_bis)^3+(-0.00239742915856412)*(T_mean_bis)*(P_mean_bis)^4+(0.000906548610701915)*(P_mean_bis)^5;
            if P_mean <=20e5
                mean_T = 3.913616357235254e+02; std_T = 28.084336211908553;  mean_P = 1.033597448034575e+06; std_P = 5.682185947619408e+05; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; cp =  (1152.30779294062)+(-14.4814709978322)*(T_mean_bis)+(112.571001285793)*(P_mean_bis)+(23.9659573259695)*(T_mean_bis)^2+(-79.4747876993283)*(T_mean_bis)*(P_mean_bis)+(40.4239532849512)*(P_mean_bis)^2+(-11.5197419972579)*(T_mean_bis)^3+(51.7212011038436)*(T_mean_bis)^2*(P_mean_bis)+(-61.8664545478244)*(T_mean_bis)*(P_mean_bis)^2+(19.9504161728721)*(P_mean_bis)^3+(2.58575382310626)*(T_mean_bis)^4+(-22.2781034965913)*(T_mean_bis)^3*(P_mean_bis)+(46.9670346989156)*(T_mean_bis)^2*(P_mean_bis)^2+(-34.0441275191149)*(T_mean_bis)*(P_mean_bis)^3+(7.67974154102382)*(P_mean_bis)^4+(-0.104769746735401)*(T_mean_bis)^5+(2.11677498777526)*(T_mean_bis)^4*(P_mean_bis)+(-10.260356659119)*(T_mean_bis)^3*(P_mean_bis)^2+(15.3326770716923)*(T_mean_bis)^2*(P_mean_bis)^3+(-7.52043313579416)*(T_mean_bis)*(P_mean_bis)^4+(0.891282637961148)*(P_mean_bis)^5;           
            else
                mean_T = 4.197338429969010e+02; std_T = 6.810759884933912;  mean_P = 2.800000000000033e+06; std_P = 4.665456720724282e+05; T_mean_bis = (T_mean - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; cp = (1848.89760028398)+(-246.117147173794)*(T_mean_bis)+(602.458603503975)*(P_mean_bis)+(190.622720489797)*(T_mean_bis)^2+(-570.375818616639)*(T_mean_bis)*(P_mean_bis)+(370.757649500479)*(P_mean_bis)^2+(-119.832148919448)*(T_mean_bis)^3+(820.961380106721)*(T_mean_bis)^2*(P_mean_bis)+(-1536.53605357561)*(T_mean_bis)*(P_mean_bis)^2+(782.289140727085)*(P_mean_bis)^3+(18.807803741737)*(T_mean_bis)^4+(-247.880350314352)*(T_mean_bis)^3*(P_mean_bis)+(1060.09591118747)*(T_mean_bis)^2*(P_mean_bis)^2+(-1732.46047881676)*(T_mean_bis)*(P_mean_bis)^3+(920.483334370913)*(P_mean_bis)^4+(-0.703636921232804)*(T_mean_bis)^5+(15.1956580638076)*(T_mean_bis)^4*(P_mean_bis)+(-116.36516230683)*(T_mean_bis)^3*(P_mean_bis)^2+(395.841792057249)*(T_mean_bis)^2*(P_mean_bis)^3+(-595.199034169249)*(T_mean_bis)*(P_mean_bis)^4+(318.397658376234)*(P_mean_bis)^5; 
            end
            Pr = cp*mu/k;
            mean_T = 4.046040468982957e+02; std_T = 24.924959340389851;  mean_P = 1.848966986417469e+06; std_P = 1.023850577380752e+06; T_mean_bis = (T_wall - mean_T)/std_T; P_mean_bis = (P_mean - mean_P)/std_P; mu_wall = (1.44569298567578e-05)+(8.04897653766163e-07)*(T_mean_bis)+(1.2046133534143e-06)*(P_mean_bis)+(3.28204169912859e-07)*(T_mean_bis)^2+(-8.13776321013716e-07)*(T_mean_bis)*(P_mean_bis)+(1.02088562667795e-06)*(P_mean_bis)^2+(-1.05166329365479e-07)*(T_mean_bis)^3+(1.25394383389519e-06)*(T_mean_bis)^2*(P_mean_bis)+(-2.48827554255973e-06)*(T_mean_bis)*(P_mean_bis)^2+(1.19193622360903e-06)*(P_mean_bis)^3+(-9.33125361157745e-09)*(T_mean_bis)^4+(-1.48694274071809e-07)*(T_mean_bis)^3*(P_mean_bis)+(1.43611038860212e-06)*(T_mean_bis)^2*(P_mean_bis)^2+(-2.76907108924106e-06)*(T_mean_bis)*(P_mean_bis)^3+(1.34110616478549e-06)*(P_mean_bis)^4+(7.14653026436247e-10)*(T_mean_bis)^5+(-9.89795305435064e-09)*(T_mean_bis)^4*(P_mean_bis)+(-4.37456744704995e-08)*(T_mean_bis)^3*(P_mean_bis)^2+(4.8157748815481e-07)*(T_mean_bis)^2*(P_mean_bis)^3+(-9.42084427674813e-07)*(T_mean_bis)*(P_mean_bis)^4+(4.76457033599083e-07)*(P_mean_bis)^5;        
            mu_rat = mu/mu_wall;

%         else
%             try
%                 mu = CoolProp.PropsSI('V',        'T', T_mean, 'P', P_mean, fluid);
%                 Pr = CoolProp.PropsSI('Prandtl',  'T', T_mean, 'P', P_mean, fluid);
%                 k  = CoolProp.PropsSI('L',        'T', T_mean, 'P', P_mean, fluid);
%                 mu_wall = CoolProp.PropsSI('V',  	'T', T_wall ,  'P', P_mean, fluid);
%                 mu_rat = mu/mu_wall;
%             catch
%                 mu = refpropm('V',          'T', T_mean,    'P', P_mean/1e3, fluid);
%                 cp = refpropm('C',          'T', T_mean,    'P', P_mean/1e3, fluid);
%                 k  = refpropm('L',          'T', T_mean,    'P', P_mean/1e3, fluid);
%                 mu_wall = refpropm('V', 	'T', T_wall ,   'P', P_mean/1e3, fluid);
%                 Pr = cp*mu/k;
%                 mu_rat = mu/mu_wall;
%             end
        end
        
    else
        if (strcmp(zone, 'liq')) && (T_mean >= T_sat_mean-5e-2)
            mu = CoolProp.PropsSI('V',        'T', T_mean, 'Q', 0, fluid);
            Pr = CoolProp.PropsSI('Prandtl',  'T', T_mean, 'Q', 0, fluid);
            k  = CoolProp.PropsSI('L',        'T', T_mean, 'Q', 0, fluid);
            mu_rat = 1;
        elseif (strcmp(zone, 'vap') || strcmp(zone, 'vap_wet') || strcmp(zone, 'tp_dryout') ) && (T_mean <= T_sat_mean +5e-2)
            mu = CoolProp.PropsSI('V',        'T', T_mean, 'Q', 1, fluid);
            Pr = CoolProp.PropsSI('Prandtl',  'T', T_mean, 'Q', 1, fluid);
            k  = CoolProp.PropsSI('L',        'T', T_mean, 'Q', 1, fluid);
            mu_rat = 1;
        else
            try
                mu = CoolProp.PropsSI('V',        'T', T_mean, 'P', P_mean, fluid);
                Pr = CoolProp.PropsSI('Prandtl',  'T', T_mean, 'P', P_mean, fluid);
                k  = CoolProp.PropsSI('L',        'T', T_mean, 'P', P_mean, fluid);
                mu_wall = CoolProp.PropsSI('V',  	'T', T_wall ,  'P', P_mean, fluid);
                mu_rat = mu/mu_wall;
            catch
                mu = refpropm('V',          'T', T_mean,    'P', P_mean/1e3, fluid);
                cp = refpropm('C',          'T', T_mean,    'P', P_mean/1e3, fluid);
                k  = refpropm('L',          'T', T_mean,    'P', P_mean/1e3, fluid);
                mu_wall = refpropm('V', 	'T', T_wall ,   'P', P_mean/1e3, fluid);
                Pr = cp*mu/k;
                mu_rat = mu/mu_wall;
            end
        end
    end
    
end
