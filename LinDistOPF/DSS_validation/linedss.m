function [] = linedss(data_location)



%% Read Text File -
branchdata_str = strcat(data_location,'\branchdata.txt');
branch = dlmread(branchdata_str);
% load 'branchdata.txt';
% branch=(branchdata);

%% LineCode and Line .dss generation
for i= 1:size(branch,1)
    a = 0;
    b = 0;
    c = 0;

    ra = [];
    rb = [];
    rc = [];
    xa = [];
    xb = [];
    xc = [];

    % RA
    if (branch(i,3)~=0)
        a = 1;
        ra = [branch(i,3)];
        xa = [branch(i,9)];
    end
    %RB
    if (branch(i,6)~=0)
        b = 1;
        rb = [branch(i,4) branch(i,6)];
        xb = [branch(i,10) branch(i,12)];
    end
    %RC
    if (branch(i,8)~=0)
        c = 1;
        rc = [branch(i,5) branch(i,7) branch(i,8)];
        xc = [branch(i,11) branch(i,13) branch(i,14)];
    end

    Line_phase(i,:) = [a b c];
    ph = a+b+c;

    if ph ==3
        R_mat = strcat(' [ ', num2str(ra),' |  ','   ',num2str(rb),' |  ','   ',num2str(rc) ,']');
        X_mat = strcat(' [ ', num2str(xa),' |  ','   ',num2str(xb),' |  ','   ',num2str(xc) ,']');
        bus_ph = '.1.2.3    ';
    end

    if ph ==2
        if (a==0)
            R_mat = strcat(' [ ', num2str(rb(2)),' |  ','   ',num2str(rc(2:3)) ,']');
            X_mat = strcat(' [ ', num2str(xb(2)),' |  ','   ',num2str(xc(2:3)) ,']');
            bus_ph = '.2.3  ';
        end
        if (b==0)
            R_mat = strcat(' [ ', num2str(ra),' |  ','   ',num2str(rc(1:2:3)) ,']');
            X_mat = strcat(' [ ', num2str(xa),' |  ','   ',num2str(xc(1:2:3)) ,']');
            bus_ph = '.1.3  ';
        end
        if (c==0)
            R_mat = strcat(' [ ', num2str(ra),' |  ','   ',num2str(rb) ,']');
            X_mat = strcat(' [ ', num2str(xa),' |  ','   ',num2str(xb) ,']');
            bus_ph = '.1.2  ';
        end
    end
    if ph ==1
        if a==1
            R_mat = strcat(' [ ', num2str(ra),']');
            X_mat = strcat(' [ ', num2str(xa),']');
            bus_ph = '.1    ';
        elseif b==1
            R_mat = strcat(' [ ', num2str(rb(2)),']');
            X_mat = strcat(' [ ', num2str(xb(2)),']');
            bus_ph = '.2    ';
        else
            R_mat = strcat(' [ ', num2str(rc(3)),']');
            X_mat = strcat(' [ ', num2str(xc(3)),']');
            bus_ph = '.3    ';
        end
    end



    if i~=1
        text_file = 'MyLineCodes.txt';
        fir = fopen(text_file, 'r');

        A = fscanf(fir,'%c');
        fclose(fir);

        text_file1 = 'Line.txt';
        fir1 = fopen(text_file1, 'r');

        B = fscanf(fir1,'%c');
        fclose(fir1);
    end
    text_file = 'MyLineCodes.txt';
    fid = fopen(text_file, 'w');

    text_file1 = 'Line.txt';
    fid1 = fopen(text_file1, 'w');

    if i~=1
        fprintf(fid,'%c' ,A);
        fprintf(fid1,'%c' ,B);
    end
    str = strcat('New linecode.',num2str(i),'     nphases=',num2str(ph),'     basefreq=60     rmatrix = ', R_mat, '               xmatrix = ', X_mat );
    fprintf(fid, str);
    fprintf(fid, '\r\n');

    str1 = strcat('New Line.L',num2str(i),'     Phases=',num2str(ph),'     Bus1 = ', num2str(branch(i,1)),bus_ph,'           Bus2 =', num2str(branch(i,2)),bus_ph, '     LineCode = ', num2str(i), '     Length = 1');
    fprintf(fid1, str1);
    fprintf(fid1, '\r\n');

    fclose('all');

end
end