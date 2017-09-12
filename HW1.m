% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num.

%your code should work no matter which of these lines is uncommented.
x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
%x = 3; y = '5'; %mixed

%your code goes here

if ischar(x)
    x = str2num(x);
end

if ischar(y)
    y = str2num(y);
end

sum = x + y;

%output your answer

fprintf('The sum is %d. \n', sum)
%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful.
% Even if you have access to the bioinformatics toolbox,
% do not use the builtin function randseq for this part.

N = 500; % define sequence length
possible_bases = 'ATGC';
rand_nums = randi(4,1,N);   % Get random numbers corresponding to bases
rand_seq = possible_bases(rand_nums);  % Assign bases to the numbers

%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF
% in your seqeunce rand_seq. Hint: see the function strfind.

% Find start and stop codons
starts = strfind(rand_seq, 'ATG');
ends = sort([strfind(rand_seq, 'TAA') strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAG')]);
longest = 0; % Initialize variable to hold longest ORF length

% Iterate through start codons
for codon = starts
    
    % Find lengths from start to all stops
    lengths = ends - codon;
    
    % Pull the shortest positive length - this is the ORF
    orf = lengths(find(lengths > 0, 1));
    
    % Save this length if it is longer than previous ones
    if orf > longest & mod(orf, 3) == 0
        % Mod 3, because ORF must act on codons 3n pase pairs away
        longest = orf;
    end
    
end

%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.

N = 500; % define sequence length
possible_bases = 'ATGC';
all_longest = zeros(1,1000);  % Initialize vector to hold all longest ORFs

% Iterate through 1000 cases
for iter = 1:1000
    rand_nums = randi(4,1,N);   % Get random numbers corresponding to bases
    rand_seq = possible_bases(rand_nums);  % Assign bases to the numbers
    
    % Find start and stop codons
    starts = strfind(rand_seq, 'ATG');
    ends = sort([strfind(rand_seq, 'TAA') strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAG')]);
    longest = 0; % Initialize variable to hold longest ORF length
    
    % Iterate through start codons
    for codon = starts
        
        % Find lengths from start to all stops
        lengths = ends - codon;
        
        % Pull the shortest positive length - this is the ORF
        orf = lengths(find(lengths > 0, 1));
        
        % Save this length if it is longer than previous ones
        if orf > longest & mod(orf, 3) == 0
            % Mod 3, because ORF must act on codons 3n pase pairs away
            longest = orf;
        end
        
    end

    % Save longest ORF lengths
    all_longest(iter) = longest;
    
end

% Determine the probability of having over 50 b.p. in the max ORF
prob_50 = length(find(all_longest > 50)) / 10;

%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 

% Define sequence lengths
seq_lens = 0:500:12000;
possible_bases = 'ATGC';
all_longest = zeros(1,1000);  % Initialize vector to hold all longest ORFs
prob_50s = zeros(1, length(seq_lens));  % Initialize vector to hold probabilities

counter = 0;
% Iterate through number of base pairs N
for N = seq_lens
    
    counter = counter + 1;
    
    % Iterate through 1000 cases
    for iter = 1:1000
        rand_nums = randi(4,1,N);   % Get random numbers corresponding to bases
        rand_seq = possible_bases(rand_nums);  % Assign bases to the numbers
        
        % Find start and stop codons
        starts = strfind(rand_seq, 'ATG');
        ends = sort([strfind(rand_seq, 'TAA') strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAG')]);
        longest = 0; % Initialize variable to hold longest ORF length
        
        % Iterate through start codons
        for codon = starts
            
            % Find lengths from start to all stops
            lengths = ends - codon;
            
            % Pull the shortest positive length - this is the ORF
            orf = lengths(find(lengths > 0, 1));
            
            % Save this length if it is longer than previous ones
            if orf > longest & mod(orf, 3) == 0
                % Mod 3, because ORF must act on codons 3n pase pairs away
                longest = orf;
            end
            
        end
        
        % Save longest ORF lengths
        all_longest(iter) = longest;
        
    end
    
    % Determine the probability of having over 50 b.p. in the max ORF
    prob_50 = length(find(all_longest > 50)) / 10;
    
    % Save probabilities
    prob_50s(counter) = prob_50;
    
end

% Plot results
plot(seq_lens, prob_50s)
ylabel('Probability of 50 b.p. ORF (% Chance)')
xlabel('Length of Sequence N')
title('Likelihood of Obtaining 50+ b.p. ORF by Sequence Length')

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

% When N is small, this curve should be near 0. This is because there is
% not much room in the sequence to even contain an ORF of length 50, and
% the codons would need to be near the start and end of the sequence.
% 
% However, when N is large, the curve should be near 100. This is because
% there are many more base pairs with which to form ORFs, so there is a
% much higher chance of getting one that is of length 50 at least.
%
% In between low and high N values, the curve should rise sharply as the
% number of base pairs increases and the possible combinations to reach 50
% base pairs rises rapidly. Then, the curve should level off and approach a
% value of 100, as the highest probability possible is 100%.
%
% The curve generated by the above code satisfies all of these criteria.

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

% Open the qPCR data
fid = fopen('qPCRdata.txt','r');

% Skip the first two lines
fgetl(fid);
fgetl(fid);

% Initialize vector to hold Cp data
Cp_data = [];

% Iterate through each line of data
for i = 1:72
    
    % Read next line and save appropriate data for Cp
    line = fgetl(fid);
    Cp_data = [Cp_data str2num(line(end-8:end-4))];
    
end

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc.

% Initialize the 6x12 array representing the plate
plate = zeros(6,12);

% Build the plate matrix from the Cp_data array
plate(1,:) = Cp_data(1:12);
plate(2,:) = Cp_data(13:24);
plate(3,:) = Cp_data(25:36);
plate(4,:) = Cp_data(37:48);
plate(5,:) = Cp_data(49:60);
plate(6,:) = Cp_data(61:72);

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

% Initialize array to hold normalized expressions
normals = zeros(6,3);

% Iterate through each condition
for condition = 1:6
    
    % Use given formula for each of the three genes
    normals(condition, 1) = 2^(mean(plate(1, 1:3)) - mean(plate(condition, 1:3)) - (mean(plate(1, 10:12)) - mean(plate(condition, 10:12))));
    normals(condition, 2) = 2^(mean(plate(1, 4:6)) - mean(plate(condition, 4:6)) - (mean(plate(1, 10:12)) - mean(plate(condition, 10:12))));
    normals(condition, 3) = 2^(mean(plate(1, 7:9)) - mean(plate(condition, 7:9)) - (mean(plate(1, 10:12)) - mean(plate(condition, 10:12))));
    
end

% Plot results
figure
labels = {'Condition 1, Gene 1' 'Condition 1, Gene 2' 'Condition 1, Gene 3' 'Condition 2, Gene 1' 'Condition 2, Gene 2' 'Condition 2, Gene 3' 'Condition 3, Gene 1' 'Condition 3, Gene 2' 'Condition 3, Gene 3' 'Condition 4, Gene 1' 'Condition 4, Gene 2' 'Condition 4, Gene 3' 'Condition 5, Gene 1' 'Condition 5, Gene 2' 'Condition 5, Gene 3' 'Condition 6, Gene 1' 'Condition 6, Gene 2' 'Condition 6, Gene 3'};
plot(1:18, [normals(1,:) normals(2,:) normals(3,:) normals(4,:) normals(5,:) normals(6,:)], 'o', 'MarkerSize', 4)
hold on
plot(1:18, ones(1,18), 'r--')
set(gca, 'XTick',1:18, 'XTickLabel',labels)
set(gca,'XTickLabelRotation',90)
xlabel('Condition and Gene')
ylabel('Normalized Gene Expression')
title('Normalized Gene Expression of qPCR Data')

%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

N = 500; % define sequence length
possible_bases = 'ATGC';
rand_nums = randi(4,1,N);   % Get random numbers corresponding to bases
rand_seq = possible_bases(rand_nums);  % Assign bases to the numbers

% Find start and stop codons
starts = strfind(rand_seq, 'ATG');
ends = sort([strfind(rand_seq, 'TAA') strfind(rand_seq, 'TGA') strfind(rand_seq, 'TAG')]);

% Get distances between all start and stop codons
diffs = bsxfun(@minus, ends, starts');

% Get indeces of possible orfs
possible_orfs = diffs > 0 & mod(diffs, 3) == 0;

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


