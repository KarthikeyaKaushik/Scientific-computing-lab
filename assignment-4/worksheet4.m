clc;
clear all;
close all;
Nxs = [7,15,31,63];
Nys = [7,15,31,63];

direct_results = zeros(2,4);
sparse_results = zeros(2,4);
gs_results = zeros(2,4);

for k = 1:4
    Nx = Nxs(k);
    Ny = Nys(k);
    
    [A,b,soln,t] = direct_solver(Nx,Ny);
    direct_results(1,k) = t;
    storage_A = whos('A');
    storage_b = whos('b');
    storage_soln = whos('soln');
    direct_results(2,k) = storage_A.bytes + storage_b.bytes + storage_soln.bytes;
    figure(k)
    subplot(211)
    surf(soln);
    title(strcat('Exact solution for : ',num2str(Nx)))
    hold on;
    figure(k+4)
    subplot(211)
    contour(soln);
    title(strcat('Exact solution for : ',num2str(Nx)))
    hold on;
    clearvars -except direct_results sparse_results gs_results k Nx Ny Nxs Nys
    
    % For sparse matrix solver
    [Asparse,b,soln_sparse,t] = sparse_solver(Nx,Ny);
    sparse_results(1,k) = t;
    storage_A = whos('Asparse');
    storage_b = whos('b');
    storage_soln = whos('soln_sparse');
    sparse_results(2,k) = storage_A.bytes + storage_b.bytes + storage_soln.bytes;
    clearvars -except direct_results sparse_results gs_results k Nx Ny Nxs Nys
    
    
    % For Gauss Seidel solver :
    [b,T,t] = gauss_seidel(Nx,Ny);
    gs_results(1,k) = t;
    storage_b = whos('b');
    storage_soln = whos('T');
    gs_results(2,k) = storage_b.bytes + storage_soln.bytes;
    figure(k)
    subplot(212)
    surf(T);
    title(strcat('Gauss Seidel solution for : ',num2str(Nx)))
    hold on;
    figure(k+4)
    subplot(212)
    contour(T);
    title(strcat('Gauss Seidel solution for : ',num2str(Nx)))

    clearvars -except direct_results sparse_results gs_results k Nx Ny Nxs Nys
    
end

T_direct = array2table(direct_results);
T_direct.Properties.RowNames = {'Run time','Storage'};
T_direct.Properties.VariableNames = {'N_7','N_15','N_31','N_63'};
T_sparse = array2table(sparse_results);
T_sparse.Properties.RowNames = {'Run time','Storage'};
T_sparse.Properties.VariableNames = {'N_7','N_15','N_31','N_63'};
T_gs = array2table(gs_results);
T_gs.Properties.RowNames = {'Run time','Storage'};
T_gs.Properties.VariableNames = {'N_7','N_15','N_31','N_63'};
disp(T_direct);
disp(T_sparse);
disp(T_gs);

%% for Gauss Seidel errors

Nxs = [7,15,31,63,127];
Nys = [7,15,31,63,127];

error_table = zeros(2,5);

for i = 1:5
    Nx = Nxs(i);
    Ny = Nys(i);
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    [~,T,~] = gauss_seidel(Nx,Ny);
    Texact = zeros(Ny,Nx);
    for j = 1:Nx
        for k = 1:Ny
            Texact(k,j) = sin(pi*k*hx)*sin(pi*j*hy);
        end
    end
    error = sqrt((sum(sum((Texact - T(2:Nx+1,2:Ny+1)).^2)))/(Nx*Ny));%error_calc(Nx,Ny,Texact,T)
    error_table(1,i) = error;
    if i == 1
        continue;
    else
        error_table(2,i) = error_table(1,i-1)/error_table(1,i);
    end
end

T_error = array2table(error_table);
T_error.Properties.RowNames = {'Error','Error_reduced'};
T_error.Properties.VariableNames = {'N_7','N_15','N_31','N_63','N_127'};
disp(T_error);




