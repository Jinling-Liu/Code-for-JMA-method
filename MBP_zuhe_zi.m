function ZUHE = MBP_zuhe_zi(zli,k)
%for the i-th component Z_l(i) of Z_l, we divide Z_l(i) to Y^1,...,Y^d, the
%possible values are donated as ZUHE
for i1=1:k
    seat_number(k-i1+1) = zli;
end
if k/2-fix(k/2)>0
    group_number=fix(k/2)+1;
else
    group_number=fix(k/2);
end

for i2=1:group_number
    if i2<group_number
        group_number_set{i2}=seat_number(1:2);
        seat_number(1:2)=[];
    else
        group_number_set{i2}=seat_number(1:end);
    end
end


for i3=1:group_number
    group_number_inset=group_number_set{i3};
    if length(group_number_set{i3})==1
        new_zuhe=[0:group_number_inset]';
    else 
        zuhe1=repmat([0:group_number_inset(1)]',group_number_inset(2)+1,1);
        re_number=repmat([0:group_number_inset(2)],group_number_inset(1)+1,1);
        total_number=(group_number_inset(2)+1)*(group_number_inset(1)+1);
        zuhe2=reshape(re_number,total_number,1);
        new_zuhe=[zuhe2 zuhe1]; new_zuhe_cj=new_zuhe*[k-i3*2+1:k-i3*2+2]';
    end
    [h_nz,~]=find(sum(new_zuhe,2)>zli+1);
    new_zuhe(h_nz,:)=[];
    zuhe_group{i3}=new_zuhe;
end

if group_number==1
    ZUHE=zuhe_group{1};
    [h_ZUHE2,~]=find(sum(ZUHE,2)==zli);
    ZUHE=ZUHE(h_ZUHE2,:);
elseif group_number==2
    [h_g2,~]=size(zuhe_group{2});
    zuhe1=repmat(zuhe_group{1},h_g2,1);
    [h_g1,~]=size(zuhe_group{1});
    zuhe2_label=repmat([1:length(zuhe_group{2})],h_g1,1);
    zuhe2_label=reshape(zuhe2_label,h_g1*h_g2,1);
    zuhe2=zuhe_group{2};
    zuhe2=zuhe2(zuhe2_label,:);
    ZUHE=[zuhe2 zuhe1];
    [h_ZUHE2,~]=find(sum(ZUHE,2)==zli);
    ZUHE=ZUHE(h_ZUHE2,:);
else
    ZUHE=zuhe_group{1};
    for i4=1:group_number-1
        [h_g2,~]=size(zuhe_group{i4+1});
        zuhe1=repmat(ZUHE,h_g2,1);
        [h_g1,~]=size(ZUHE);
        zuhe2_label=repmat([1:h_g2],h_g1,1);
        zuhe2_label=reshape(zuhe2_label,h_g2*h_g1,1);
        zuhe2=zuhe_group{i4+1};
        zuhe2=zuhe2(zuhe2_label,:);
        ZUHE=[zuhe2 zuhe1];
        [h_ZUHE2,~]=find(sum(ZUHE,2)<zli+1);
        ZUHE=ZUHE(h_ZUHE2,:);
    end
    [h_ZUHE2,~]=find(sum(ZUHE,2)==zli);
    ZUHE=ZUHE(h_ZUHE2,:);
end
