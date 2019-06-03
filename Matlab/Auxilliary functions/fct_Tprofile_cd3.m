function out = fct_Tprofile_cd3(dt_Tol_tp, dT_tol_grad,IRP_path, filename, photo_IR1, photo_IR2, photo_IR3, photo_IR4, T_sat_cd_su_ref, T_sat_cd_ex_ref, T_cd_su, T_cd_contact_fh, T_cd_contact_bh,  T_cd_contact_fl,  T_cd_contact_bl,  T_cd_ex, lim_liq_man, lim_vap_man)

out.flag = 1;

%% POSITION DEFINITION
out.vec_pos_ref = 0:14;


%% CHECK PHOTO RESSOURCE
out.take_photo_1 = min(1,exist([IRP_path 'IRP_' filename '\IRP_' num2str(photo_IR1) '_results.mat'],'file'));
out.take_photo_2 = min(1,exist([IRP_path 'IRP_' filename '\IRP_' num2str(photo_IR2) '_results.mat'],'file'));
out.take_photo_3 = min(1,exist([IRP_path 'IRP_' filename '\IRP_' num2str(photo_IR3) '_results.mat'],'file'));
out.take_photo_4 = min(1,exist([IRP_path 'IRP_' filename '\IRP_' num2str(photo_IR4) '_results.mat'],'file'));

%% SATURATION CONDITIONS
out.T_sat_cd_su_ref = T_sat_cd_su_ref;
out.T_sat_cd_ex_ref = T_sat_cd_ex_ref;

%% CONTACT THERMOCOUPLES
out.T_cd_contact(1,:) = [T_cd_su T_cd_contact_fh T_cd_contact_bh  T_cd_contact_fl  T_cd_contact_bl  T_cd_ex];
out.vec_pos_contact = [-1 2 5 8 11 15];

%% IR TEMPERATURE RECORDING Photo 2 - zoom front view
if out.take_photo_2
    load([IRP_path 'IRP_' filename '\IRP_' num2str(photo_IR2) '_results.mat'], ['IRP_' num2str(photo_IR2) '_num'],['IRP_' num2str(photo_IR2) '_tube_2_PosVec'],['IRP_' num2str(photo_IR2) '_tube_5_PosVec'])
    eval(['out.IR_temp_photo2 = IRP_' num2str(photo_IR2) '_num;']);
    eval(['IR_pos_t2_photo_2_bis = IRP_' num2str(photo_IR2) '_tube_2_PosVec;']);
    eval(['IR_pos_t5_photo_2_bis = IRP_' num2str(photo_IR2) '_tube_5_PosVec;']);
    out.IR_pos_t5_photo_2(:,1) = fliplr(IR_pos_t5_photo_2_bis(:,1)')'; out.IR_pos_t5_photo_2(:,2) = fliplr(IR_pos_t5_photo_2_bis(:,2)')';
    out.IR_pos_t2_photo_2(:,1) = fliplr(IR_pos_t2_photo_2_bis(:,1)')'; out.IR_pos_t2_photo_2(:,2) = fliplr(IR_pos_t2_photo_2_bis(:,2)')';
    for i_t2 = 1:length(out.IR_pos_t2_photo_2(:,1))
        if sum(isnan(out.IR_pos_t2_photo_2(i_t2,:)))
            out.T_ir_front_t2(1,i_t2) = NaN;
        else
            out.T_ir_front_t2(1,i_t2) =  out.IR_temp_photo2(out.IR_pos_t2_photo_2(i_t2,2),out.IR_pos_t2_photo_2(i_t2,1));
        end
    end
    for i_t5 = 1:length(out.IR_pos_t5_photo_2(:,1))
        if sum(isnan(out.IR_pos_t5_photo_2(i_t5,:)))
            out.T_ir_front_t5(1,i_t5) = NaN;
        else
            out.T_ir_front_t5(1,i_t5) =  out.IR_temp_photo2(out.IR_pos_t5_photo_2(i_t5,2),out.IR_pos_t5_photo_2(i_t5,1));
        end
    end
    
else
    out.T_ir_front_t5 = NaN*ones(1,6);
    out.T_ir_front_t2 = NaN*ones(1,6);
    out.IR_temp_photo2 = NaN;
    out.IR_pos_t5_photo_2 = NaN;
    out.IR_pos_t2_photo_2 = NaN;
end
out.vec_pos_front = [2 4 6 8 10 12];


%% IR TEMPERATURE RECORDING Photo 3 - zoom back  view
if out.take_photo_3
    load([IRP_path 'IRP_' filename '\IRP_' num2str(photo_IR3) '_results.mat'],['IRP_' num2str(photo_IR3) '_num'],['IRP_' num2str(photo_IR3) '_tube_2_PosVec'],['IRP_' num2str(photo_IR3) '_tube_3_PosVec'])
    eval(['out.IR_temp_photo3 = IRP_' num2str(photo_IR3) '_num;']);
    eval(['IR_pos_t2_photo_3_bis = IRP_' num2str(photo_IR3) '_tube_2_PosVec;']);
    eval(['IR_pos_t3_photo_3_bis = IRP_' num2str(photo_IR3) '_tube_3_PosVec;']);
    out.IR_pos_t2_photo_3(:,1) = fliplr(IR_pos_t2_photo_3_bis(:,1)')'; out.IR_pos_t2_photo_3(:,2) = fliplr(IR_pos_t2_photo_3_bis(:,2)')';
    out.IR_pos_t3_photo_3(:,1) = fliplr(IR_pos_t3_photo_3_bis(:,1)')'; out.IR_pos_t3_photo_3(:,2) = fliplr(IR_pos_t3_photo_3_bis(:,2)')';
    
    for i_t2 = 1:length(out.IR_pos_t2_photo_3(:,1))
        if sum(isnan(out.IR_pos_t2_photo_3(i_t2,:)))
            out.T_ir_back_t2(1,i_t2) = NaN;
        else
            out.T_ir_back_t2(1,i_t2) =  out.IR_temp_photo3(out.IR_pos_t2_photo_3(i_t2,2),out.IR_pos_t2_photo_3(i_t2,1));
        end
    end
    for i_t3 = 1:length(out.IR_pos_t3_photo_3(:,1))
        if sum(isnan(out.IR_pos_t3_photo_3(i_t3,:)))
            out.T_ir_back_t3(1,i_t3) = NaN;
        else
            out.T_ir_back_t3(1,i_t3) =  out.IR_temp_photo3(out.IR_pos_t3_photo_3(i_t3,2),out.IR_pos_t3_photo_3(i_t3,1));
        end
    end
else
    out.T_ir_back_t3 = NaN*ones(1,9);
    out.T_ir_back_t2 = NaN*ones(1,9);
    out.IR_temp_photo3 = NaN;
    out.IR_pos_t3_photo_3 = NaN;
    out.IR_pos_t2_photo_3 = NaN;
end
out.vec_pos_back = [0   1   3   5   7   9   11  13   14];

%% IR TEMPERATURE RECORDING Photo 4 - side back  view
if out.take_photo_4
    load([IRP_path 'IRP_' filename '\IRP_' num2str(photo_IR4) '_results.mat'],['IRP_' num2str(photo_IR4) '_num'],['IRP_' num2str(photo_IR4) '_tube_X_PosVec'])
    eval(['out.IR_temp_photo4 = IRP_' num2str(photo_IR4) '_num;']);
    eval(['out.IR_pos_photo_4 = IRP_' num2str(photo_IR4) '_tube_X_PosVec;']);
    
    for i_tx = 1:length(out.IR_pos_photo_4(:,1))
        if sum(isnan(out.IR_pos_photo_4(i_tx,:)))
            out.T_ir_back_side_tx(1,i_tx) = NaN;
        else
            out.T_ir_back_side_tx(1,i_tx) =  out.IR_temp_photo4(out.IR_pos_photo_4(i_tx,2),out.IR_pos_photo_4(i_tx,1));
        end
    end
else
    out.T_ir_back_side_tx(1,:) = NaN*ones(1,9);
    out.IR_pos_photo_4 = NaN;
    out.IR_temp_photo4 = NaN;
end



%CORRECTION 1 : Assemble Photo 4 and photo 3 -> create T_ref_back
if out.take_photo_4 && out.take_photo_3
    vec2compare_photo3_photo4 = [4 5 6 7];
    vec2compare_photo3_photo4 = vec2compare_photo3_photo4(not(isnan(out.T_ir_back_t2(1,vec2compare_photo3_photo4))) & not(isnan(out.T_ir_back_side_tx(1,vec2compare_photo3_photo4))));
    difference_photo_3_photo_4 = out.T_ir_back_side_tx(1, vec2compare_photo3_photo4) - out.T_ir_back_t2(1, vec2compare_photo3_photo4);
    shift_photo3_photo_4 = mean(difference_photo_3_photo_4);
    out.T_ir_back_side_tx_bis(1,:) = out.T_ir_back_side_tx(1,:) - shift_photo3_photo_4;
    out.T_ir_ref_back(1, [1 3 4 5 6 7 9]) = out.T_ir_back_t2(1, [1 3 4 5 6 7 9]);
    if isnan(out.T_ir_back_side_tx_bis(1, 2))
        out.T_ir_ref_back(1, 2) = out.T_ir_back_t2(1, 2);
    else
        out.T_ir_ref_back(1, 2) = out.T_ir_back_side_tx_bis(1, 2);
    end
    if isnan(out.T_ir_back_side_tx_bis(1, 8))
        out.T_ir_ref_back(1, 8) = out.T_ir_back_t2(1, 8);
    else
        out.T_ir_ref_back(1, 8) = out.T_ir_back_side_tx_bis(1, 8);
    end
    
else
    out.T_ir_ref_back(1,:) = out.T_ir_back_t2(1, :);
    out.T_ir_back_side_tx_bis = out.T_ir_back_side_tx;
end

%CORRECTION 2 : fisrt assembly of front and back  to ensure linearity in the middle of the condenser
if out.take_photo_2 && out.take_photo_3
    pos2compare_front_back = [7 9];
    T_ir_front_fit = fit(out.vec_pos_front',out.T_ir_front_t2(1,:)','linear');
    shift_front_back = mean(T_ir_front_fit(pos2compare_front_back)' - out.T_ir_ref_back(1,ismember(out.vec_pos_back, pos2compare_front_back)));
    out.T_ir_ref_front(1, :) = out.T_ir_front_t2(1, :) - shift_front_back ;
    out.T_ir_ref(1, [1 2 4 6 8 10 12 14 15]) =  out.T_ir_ref_back(1,:);
    out.T_ir_ref(1, [3 5 7 9 11 13]) =  out.T_ir_ref_front(1,:);
else
    out.flag = -2;
    out.T_ir_ref_front = out.T_ir_front_t2;
    out.T_ir_ref(1, :) = NaN*ones(1,15);
end

%CORRECTION 2 : final correction to have proper location of the two-phase region on the both sides
if out.take_photo_2 && out.take_photo_3
    vec2analyse_tp = 1:14;
    diff_Tsat_Tir_ref = abs(out.T_ir_ref(1,vec2analyse_tp)-out.T_sat_cd_su_ref);
    i_small_dev_Tsat = vec2analyse_tp((diff_Tsat_Tir_ref<5));
    
    diff_Tinter_ir_ref = abs(diff(out.T_ir_ref(1,i_small_dev_Tsat)));
    i_small_dev_inter = i_small_dev_Tsat((diff_Tinter_ir_ref<1.5));
    if length(i_small_dev_inter) > 1
        i_small_dev_inter = i_small_dev_inter(not(ismember(i_small_dev_inter, [1])));
    end
    Nbr_zone_tp = length(i_small_dev_inter)+1;
    postion_start_tp = i_small_dev_inter(1);
    T_start_tp(1,1) = out.T_ir_ref(1,postion_start_tp);
    if ismember(postion_start_tp, [1 2 4 6 8 10 12 14 15])
        side_start_tp = 'back';
        shift_back_2 = out.T_ir_ref(1,postion_start_tp)-out.T_sat_cd_su_ref;
        out.T_ir_ref_2(1, [1 2 4 6 8 10 12 14 15]) = out.T_ir_ref(1, [1 2 4 6 8 10 12 14 15]) - shift_back_2;
        if Nbr_zone_tp>2
            shift_front_2 = out.T_ir_ref(1,postion_start_tp+1) - (0.5*out.T_ir_ref_2(1,postion_start_tp)+0.5*out.T_ir_ref_2(1,postion_start_tp+2));
            
        else
            shift_front_2 = out.T_ir_ref(1,postion_start_tp+1) - out.T_ir_ref_2(1,postion_start_tp);
        end
        out.T_ir_ref_2(1, [3 5 7 9 11 13]) = out.T_ir_ref(1, [3 5 7 9 11 13]) - shift_front_2;
    elseif ismember(postion_start_tp, [3 5 7 9 11 13])
        side_start_tp = 'front';
        shift_front_2 = out.T_ir_ref(1,postion_start_tp)-out.T_sat_cd_su_ref;
        out.T_ir_ref_2(1, [3 5 7 9 11 13]) = out.T_ir_ref(1, [3 5 7 9 11 13]) - shift_front_2;
        if Nbr_zone_tp>2
            shift_back_2 = out.T_ir_ref(1,postion_start_tp+1) - (0.5*out.T_ir_ref_2(1,postion_start_tp)+0.5*out.T_ir_ref_2(1,postion_start_tp+2));
            
        else
            shift_back_2 = out.T_ir_ref(1,postion_start_tp+1) - out.T_ir_ref_2(1,postion_start_tp);
        end
        out.T_ir_ref_2(1, [1 2 4 6 8 10 12 14 15]) = out.T_ir_ref(1, [1 2 4 6 8 10 12 14 15]) - shift_back_2;
    end
    out.final_shif_front = out.T_ir_ref_2(1,3) - out.T_ir_front_t2(1, 1);
    out.final_shif_back  = out.T_ir_ref_2(1,6) - out.T_ir_back_t2(1, 4);
else
    out.final_shif_front = NaN;
    out.final_shif_back = NaN;
    out.T_ir_ref_2(1, :) = NaN*ones(1,15);
end

out.T_ir_ref_2(1, [1 2]) =      [T_cd_su T_cd_su];
out.T_ir_ref_2(1, [14 15]) =    [T_cd_ex T_cd_ex];

%% ZONE DEFINITION
if out.take_photo_2 && out.take_photo_3
    
    Nbr_zone = length(out.T_ir_ref_2(1,:))-1;
    
    %Criteria based on proximity of Tsat
    for i_zone = 1:Nbr_zone
        if ((out.T_ir_ref_2(1,i_zone) + out.T_ir_ref_2(1,i_zone+1))/2) < out.T_sat_cd_su_ref + dt_Tol_tp && ((out.T_ir_ref_2(1,i_zone) + out.T_ir_ref_2(1,i_zone+1))/2) > out.T_sat_cd_ex_ref - dt_Tol_tp
            crit_tp1(i_zone) = 1 ;
        else
            crit_tp1(i_zone) = 0 ;
        end
    end
    
    % Criteria based on gradient of temperature
    diff_T_ir_ref_2 = diff(out.T_ir_ref_2(1,:));
    for i_zone = 1:Nbr_zone
        if diff_T_ir_ref_2(i_zone) < dT_tol_grad
            crit_tp2(i_zone) = 0 ;
        else
            crit_tp2(i_zone) = 1 ;
        end
    end
    
    
    
    if isempty(lim_vap_man) || isempty(lim_liq_man)
        zone_tp = find(crit_tp1==1 & crit_tp2==1);
        
        %Criteria of continuity of the two phase region continuity
        zone_tp = zone_tp(1):zone_tp(end);
        
        % Criteria of logic zone sequence
        zone_vap = 1:zone_tp(1)-1;
        zone_liq = zone_tp(end)+1:Nbr_zone;
        if not(isempty(zone_tp))
            for i_tp = zone_tp
                out.zone_cd{1,i_tp} = 'tp';
            end
        end
        if not(isempty(zone_vap))
            for i_vap = zone_vap
                out.zone_cd{1,i_vap} = 'vap';
            end
        end
        if not(isempty(zone_liq))
            for i_liq = zone_liq
                out.zone_cd{1,i_liq} = 'liq';
            end
        end
    else
        
        % Manual setting of zones based on external input
        zone_tp_man = lim_vap_man:lim_liq_man;
        zone_vap_man = 1:zone_tp_man(1)-1;
        zone_liq_man = zone_tp_man(end)+1:Nbr_zone;
        if not(isempty(zone_tp_man))
            for i_tp = zone_tp_man
                out.zone_cd{1,i_tp} = 'tp';
            end
        end
        if not(isempty(zone_vap_man))
            for i_vap = zone_vap_man
                out.zone_cd{1,i_vap} = 'vap';
            end
        end
        if not(isempty(zone_liq_man))
            for i_liq = zone_liq_man
                out.zone_cd{1,i_liq} = 'liq';
            end
        end
    end
    
else
    out.zone_cd = '';
end

out = orderfields(out);
end