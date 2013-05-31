clear all, close all, clc

% load files from 0 -> mpi.size

maxSize = 128;
numOpArray = zeros(maxSize,maxSize);


filename = 'numOp_rank';
filetype = 'txt';

% loop files
for i = 0 : maxSize-1
    % read formatted file
    fid = fopen(['stout/',filename,num2str(i),'_830268.txt']);
    if fid==-1 % file not found
    else
        data = textscan(fid, '%d %d %u64'); % myrank, mpi.size, myCounter
        fclose(fid);
        numOpArray(i +1,data{2}(:)) = data{3}(:)';
        numSizeArray(i +1,data{2}(:)) = data{2}(:)';
    end
end

%% sum
totOpVec = sum(numOpArray,1);

%% plot
figure(1),clf
[I,J] = find(numSizeArray);
plot(numSizeArray(I,J),numOpArray(I,J),'x')
xlabel('mpi.size')
ylabel('"number of loops" per process')
grid on
set(gca,'XScale','log','YScale','log')

figure(2),clf
plot(numSizeArray(1,numSizeArray(1,:)>0),totOpVec(numSizeArray(1,:)>0),'x')
xlabel('mpi.size')
ylabel('total "number of loops"')
grid on
set(gca,'XScale','lin','YScale','log')