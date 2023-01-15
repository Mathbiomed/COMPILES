% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    steadyState                                                            %
%                                                                           %
%                                                                           %
% OUTPUT: Returns the steady state solution of a chemical reaction network  %
%    (CRN) parametrized by rate constants (and/or sigma's). The free        %
%    parameters and conservation laws are also listed after the solution.   %
%    If there are subnetworks that could not be solved because translation  %
%    is taking too long, solving the subnetwork is skipped and a message    %
%    saying it could not be solved is displayed. In this case, the          %
%    parametrization of the steady state of the entire network is not       %
%    completed but the lists of solved and unsolved subnetworks are         %
%    displayed. In the case where all subnetworks are solved but the        %
%    solution of the entire network cannot be parametrized in terms of the  %
%    free parameters (due to two or more species dependent on each other),  %
%    the solution returned is in terms of both free parameters and other    %
%    "nonfree" species. If there are subnetworks containing only 1          %
%    reaction, it means that there are species with 0 steady steady; hence, %
%    the network has not positive steady state and a message appears saying %
%    so. The output variables 'equation', 'species', 'free_parameter',      %
%    'conservation_law', and 'model' allow the user to view the following,  %
%    respectively:                                                          %
%       - List of parametrization of the steady state of the system         %
%       - List of steady state species of the network                       %
%       - List of free parameters of the steady state                       %
%       - List of conservation laws of the system                           %
%       - Complete network with all the species listed in the 'species'     %
%            field of the structure 'model'                                 %
%                                                                           %
% INPUT: model: a structure, representing the CRN (see README.txt for       %
%    details on how to fill out the structure)                              %
%                                                                           %
% Notes:                                                                    %
%    1. It is assumed that the CRN has mass action kinetics.                %
%    2. The algorithm first decomposes the network intro its finest         %
%          independent decomposition. Then the steady state of the species  %
%          in each subnetwork is solved. Finally, these solutions are       %
%          combined to get the solution to the entire network.              %
%    3. Decomposition into the finest independent decomposition comes from  %
%          [2].                                                             %
%    4. Translation of network comes from [3].                              %
%    5. Parametrization of steady state solution comes from [4].            %
%    6. Computation of directed spanning trees toward vertices come from    %
%          [1].                                                             %
%    7. Sometimes, when the solution could not be expressed in terms of     %
%          free parameters only, renaming the variables can solve the       %
%          problem. Some subnetworks may be solved for different species    %
%          depending on the variable assigned to them since the selection   %
%          of species to solve is based on alphabetical order.              %
%    8. Ideas for some parts of the code was motivated by [5].              %
%                                                                           %
% References                                                                %
%    [1] Gabow H and Meyers E (1978) Finding all spanning trees of directed %
%           and undirected graphs. SIAM J Comput 7(3):280-287.              %
%           https://doi.org/10.1137/0207024                                 %
%    [2] Hernandez B, De la Cruz R (2021) Independent decompositions of     %
%           chemical reaction networks. Bull Math Biol 83(76):1–23.         %
%           https://doi.org/10.1007/s11538-021-00906-3                      %
%    [3] Hong H, Hernandez B, Kim J, Kim JK (2022) Computational            %
%           translation framework identifies biochemical reaction networks  %
%           with special topologies and their long-term dynamics            %
%           (submitted)                                                     %
%    [4] Johnston M, Mueller S, Pantea C (2019) A deficiency-based approach %
%           to parametrizing positive equilibria of biochemical reaction    %
%           systems. Bull Math Biol 81:1143–1172.                           %
%           https://doi.org/10.1007/s11538-018-00562-0                      %
%    [5] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical        %
%           reaction network theory. Bioinform 25(21):2853–2854.            %
%           https://doi.org/10.1093/bioinformatics/btp513                   %
%                                                                           %
% Created: 15 July 2022                                                     %
% Last Modified: 28 September 2022                                          %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [equation, species, free_parameter, conservation_law, model] = steadyState(model)

    %
    % Step 1: Decompose the network into its finest independent decomposition
    %

    [model, N, R, G, P] = indepDecomp(model);
    
    if numel(P) == 1
        fprintf('The network has no nontrival independent decomposition.\n\n')
    else
        fprintf('The network has %d subnetworks.\n\n', numel(P))
    end
    
    % Create a vector of reaction numbers for the total number of reactions
    reac_num = [ ];
    for i = 1:numel(model.reaction)
        if model.reaction(i).reversible == 0
            reac_num(end+1) = i;
        else
            reac_num(end+1) = i;
            reac_num(end+1) = i;
        end
    end

    % Get the number of reactions in each subnetwork
    subnetwork_reaction_count = cellfun(@numel, P);

    % Look for subnetworks with just 1 reaction
    single_reaction = find(subnetwork_reaction_count == 1);

    % Exit the algorithm if there is a subnetwork with just 1 reaction
    if ~isempty(single_reaction)
        if numel(single_reaction) == 1
            equation = { };
            species = [ ];
            free_parameter = [ ];
            conservation_law = { };
            fprintf('Subnetwork %s has 1 reaction only.\n', single_reaction)
            fprintf('The network has no positive steady state.\n\n')
        else
            equation = { };
            species = [ ];
            free_parameter = [ ];
            conservation_law = { };
            print_list = sprintf('%d, ', single_reaction);
            print_list(end-1:end) = [ ]; % To clean the trailing comma and space at the end of the list
            fprintf('Subnetworks %s have 1 reaction only.\n', print_list)
            fprintf('The network has no positive steady state.\n\n')
        end
        return
    end



    %
    % Step 2: Get the conservation laws of the system of ODEs
    %

    % Get the conservation laws: the kernel of the transpose of the stoichiometric subspace
    conservation_law_matrix = null(N', 'r');


    % Initialize list of conservation laws
    conservation_law = { };

    % Get the conservation laws
    for i = 1:size(conservation_law_matrix, 2)

        % Get the column
        conservation_law_ = conservation_law_matrix(:, i);

        % Locate the index of nonzero entries of each column
        conservation_law_nnz = find(conservation_law_);

        % Get the conservation law
        for j = 1:length(conservation_law_nnz)
        
            % First entry
            if j == 1

                % Do not show coefficient if it is 1
                if conservation_law_(conservation_law_nnz(j)) == 1
                    law = ['d' model.species{conservation_law_nnz(j)} '/dt'];

                % Just show a negative sign if it is -1
                elseif conservation_law_(conservation_law_nnz(j)) == -1
                    law = ['-d' model.species{conservation_law_nnz(j)} '/dt'];
                
                % All other cases
                else
                    law = [num2str(conservation_law_(conservation_law_nnz(j))) 'd' model.species{conservation_law_nnz(j)} '/dt'];
                end

            % Succeeding entries
            else

                % Do not show coefficient if it is 1
                if conservation_law_(conservation_law_nnz(j)) == 1
                    law = [law ' + d' model.species{conservation_law_nnz(j)} '/dt'];

                % For negative coefficients
                elseif conservation_law_(conservation_law_nnz(j)) < 0

                    % Just show a negative sign if it is -1
                    if conservation_law_(conservation_law_nnz(j)) == -1
                        law = [law ' - d' model.species{conservation_law_nnz(j)} '/dt'];
                    
                    % All other negative coefficients
                    else
                        law = [law ' - ' num2str(abs(conservation_law_(conservation_law_nnz(j)))) 'd' model.species{conservation_law_nnz(j)} '/dt'];
                    end

                % All other cases
                else
                    law = [law ' + ' num2str(conservation_law_(conservation_law_nnz(j))) 'd' model.species{conservation_law_nnz(j)} '/dt'];
                end
            end
        end

        % Add the law to the list
        conservation_law{end+1} = [law ' = 0'];
    end



    %
    % Step 3: Solve for the steady state of each subnetwork
    %
    
    % Initialize list of constants used in parametrization
    B_ = [ ];
    sigma_ = [ ];
    
    % Initialize list of equations of the steady state solution
    equation = { };             % List of all steady state solutions
    equation_for_solving = { }; % Same list but formatted for using the 'solve' function
    
    % Initialize list of steady state species
    species = [ ];
    
    % Initialize list of unsolved subnetworks
    unsolved_subnetwork = [ ];

    % Initialize checker if the entire network needs to be solved
    % Entire network will be solved only if all subnetworks are solved
    do_not_solve_entire_network = 0;

    % Go through each subnetwork
    for i = 1:numel(P)

        % Get the corresponding reactions of the subnetwork from the parent network
        model_P(i).id = [model.id ' - Subnetwork ' num2str(i)];
        model_P(i).species = { };
        reac_P = unique(reac_num(P{i}));
        for j = 1:numel(reac_P)
            model_P(i).reaction(j) = model.reaction(reac_P(j));
        end

        % Create a vector of reaction numbers for the total number of reactions in the subnetwork
        reac_num_ = [ ];
        for j = 1:numel(model_P(i).reaction)
            if model_P(i).reaction(j).reversible == 0
                reac_num_(end+1) = j;
            else
                reac_num_(end+1) = j;
                reac_num_(end+1) = j;
            end
        end

        % Get the unique reaction numbers
        reac_num_unique = unique(reac_num_);

        % Initialize list of reactions of the subnetwork
        reac_subnetwork = { };

        % Create list of reactions of the subnetwork
        for j = 1:numel(reac_num_unique)

            % If the reaction is reversible
            if model_P(i).reaction(reac_num_unique(j)).reversible == 1

                % Split the reactant and product complexes
                split_complex = split(model_P(i).reaction(reac_num_unique(j)).id, '<->');
                
                % Add to the list the two reactions separately
                reac_subnetwork{end+1} = [split_complex{1}, '->', split_complex{2}];
                reac_subnetwork{end+1} = [split_complex{2}, '->', split_complex{1}];
            else
                reac_subnetwork{end+1} = model_P(i).reaction(reac_num_unique(j)).id;
            end
        end

        % Get the reaction numbers from the independent decomposition
        P_ = P{i};
    
        if numel(P) == 1
            fprintf('- Network -\n\n')
        else
            fprintf('- Subnetwork %d -\n\n', i)
        end

        % Show reactions of the subnetwork
        for j = 1:numel(P_)
            fprintf('R%d: %s\n', P_(j), reac_subnetwork{j})
        end

        if numel(P) == 1
            fprintf('\nSolving the network...\n\n')
        else
            fprintf('\nSolving Subnetwork %d...\n\n', i)
        end

        % Solve for the the steady state of the subnetwork
        [param, B_, sigma_, k, sigma, tau, model_P(i), skip] = analyticSolution(model_P(i), P_, B_, sigma_);

        % If finding a translation takes too long for the subnetwork
        if skip == 1

            % Case 1: No nontrival independent decomposition
            if numel(P) == 1
                fprintf('Could not solve the network.\n\n')
                
                % List the subnetwork as unsolved
                unsolved_subnetwork(end+1) = i;

                % No free parameters found (since we could not solve the network)
                free_parameter = [ ];

                % Let the checker know that we don't need to solve the entire network
                do_not_solve_entire_network = 1;

                % Exit the loop
                break

            % Case 2: There are other subnetworks
            else
                fprintf('Could not solve Subnetwork %d.\n\n', i)

                % List the subnetwork as unsolved
                unsolved_subnetwork(end+1) = i;

                % Let the checker know that we don't need to solve the entire network
                do_not_solve_entire_network = 1;

                % Go to the next iteration already
                continue
            end
        end

        fprintf('\n')
    
        % Substitute the tau's for the free parameters of the subnetwork
        for j = 1:length(param)
            if ismember(param(j), tau)
                param = subs(param, param(j), model_P(i).species{j});
            end
        end

        % Substitute the sigma's for the free parameters of the subnetwork
        for j = 1:length(model_P(i).species)
            
            % Get the numerator and denominator of each solution
            [N, D] = numden(param(j));
    
            % Check if it is of the form sigma/k or k/sigma
            if ismember(N, sigma) | ismember(D, sigma)
    
                % Solve for sigma in terms of k, then substitute
                if ismember(N, sigma)
                    sigma_solved = solve(param(j) == model_P(i).species{j}, sigma(find(ismember(sigma, N))));
                    param = subs(param, N, sigma_solved);
                elseif ismember(D, sigma)
                    sigma_solved = solve(param(j) == model_P(i).species{j}, sigma(find(ismember(sigma, D))));
                    param = subs(param, D, sigma_solved);
                end
            end
        end

        % Go through each species of the subnetwork
        for j = 1:length(model_P(i).species)

            % Get only the steady state species
            if model_P(i).species{j} ~= string(param(j))
                equation{end+1} = param(j);
                equation_for_solving{end+1} = model_P(i).species{j} == param(j);
                
                % Collect all steady state species of the subnetwork
                species = [species, cell2sym(model_P(i).species(j))];
            end
        end
    end



    %
    % Step 4: Use the solutions of each subnetwork to get the steady state of the entire network
    %

    % If the checker says we don't need to solve the entire network
    if do_not_solve_entire_network == 1

        % No free parameters found (since we could not solve the network)
        free_parameter = [ ];

        % Get list of solved subnetworks
        solved = setdiff(1:numel(P), unsolved_subnetwork);
        print_solved = sprintf('%d, ', solved);
        print_solved(end-1:end) = [ ]; % To clean the trailing comma and space at the end of the list
        if isempty(solved)
            if numel(P) == 1
                fprintf('Could not solve the network.\n')
            else
                fprintf('No subnetwork was solved.\n')
            end
        elseif length(solved) == 1
            fprintf('Solved only Subnetwork %s.\n', print_solved)
        else
            fprintf('Solved only Subnetworks %s.\n', print_solved)
        end

        % Get list of unsolved subnetworks
        print_unsolved = sprintf('%d, ', unsolved_subnetwork);
        print_unsolved(end-1:end) = [ ]; % To clean the trailing comma and space at the end of the list
        if isempty(solved)
            fprintf('Try to compute the subnetworks manually then combine their solutions.\n\n')
        elseif length(unsolved_subnetwork) == 1
            if numel(P) == 1
                fprintf('Try to compute the network manually.\n\n')
            else
                fprintf('Try to compute Subnetwork %s manually then combine with previously solved subnetworks.\n\n', print_unsolved)
            end
        elseif isempty(solved) & length(unsolved_subnetwork) > 1
            fprintf('Try to compute Subnetworks %s manually then combine their solutions.\n\n', print_unsolved)
        else
            fprintf('Try to compute Subnetworks %s manually then combine with previously solved subnetworks.\n\n', print_unsolved)
        end

        % Display conservation laws
        if isempty(conservation_law_matrix)
            fprintf('Conservation laws: None \n\n')
        elseif length(conservation_law) == 1
            fprintf('Conservation law: %s \n\n', conservation_law{1})
        else
            fprintf('Conservation laws: \n')
            for i = 1:length(conservation_law)
                disp(conservation_law{i})
            end
        end

        % Exit the algorithm
        return
    end

    if numel(P) == 1
        fprintf('Solving positive steady state parametrization of the network...\n\n')
    else
        fprintf('Solving positive steady state parametrization of the entire network...\n\n')
    end
    
    % Get the unique species we have solved so far
    [unique_species, ~, unique_numbering] = unique(species);
    
    % Create a list of species already solved or reserved to be solved
    species_solved = unique_species;

    % Count how many of each species appeared in the initial solution
    unique_species_occurrence = histc(unique_numbering, unique(unique_numbering));
    
    % Get index of species that occur more than once
    not_unique_species_index = find(unique_species_occurrence > 1);
    
    % Go through each species with duplicate equation
    for i = 1:length(not_unique_species_index)
    
        % Look for the equations of the species
        equation_number = find(has(species, unique_species(not_unique_species_index(i))));
        
        % Initialize list of species in each equation
        species_ = [ ];
    
        % Get all species in each equation
        for j = 1:length(equation_number)
            
            % Get all variables involved
            all_var = symvar(equation{equation_number(j)});
    
            % Get the species only
            species_ = [species_, intersect(cell2sym(model.species), all_var)];
        end

        % Removed from considerations the reserved species
        species_ = setdiff(species_, species_solved);
        
        % List down the variables to be solved: should be the number of equations
        % One of them should be the repeated species
        % The rest will come from 'species_'
        vars = [unique_species(not_unique_species_index(i))];
    
        % Determine the possible additional species to be solved for
        additional = nchoosek(species_, length(equation_number) - 1);
    
        % Solve the system of equations
        for j = 1:length(additional)
            
            % Add each combination of possible additional species
            vars_ = [vars, additional(j,:)];

            % Solve for the indicated species
            solution = solve([equation_for_solving{equation_number}], vars_);
    
            % Convert the structure of solutions to cell array
            solution = struct2cell(solution);

            % Initialize solution_
            solution_ = { };

            % For multiple solutions, we'll just get the first one
            for i = 1:length(solution)
                solution_{end+1} = solution{i}(1);
            end

            % This is the solution we'll get
            solution = solution_;

            % Check if the solution is not empty
            if ~isempty(cell2sym(solution))
    
                % Prepare a checker variable
                checker = [ ];

                % Check if any of the solutions has a negative sign
                for k = 1:length(solution)
                    checker = [checker, findstr(string(solution{k}), '-')];
                end

                % No need to consider the other combinations of possible additional species
                if isempty(checker)
                    break
                end
            end
        end

        % Add the species you just solved to the list of solved species
        species_solved = unique([species_solved, vars_]);

        % Get the species of the solution and replace these to the list of species
        species(equation_number) = vars_;

        % Get the solutions and replace these to the list of equations
        equation(equation_number) = solution;
    end

    % Initialize list of equations in terms of the free parameters
    equation_final = equation;

    % Initialize flag whether to skip the simplification of the solution
    skip_simplification = 0;
    
    % Make sure the equations contain only free parameters
    for i = 1:length(equation_final)

        % Initialize list of dependent variables
        dep_var_cycle = {species(i)};
    
        % Keep on substituting until there are no more steady state variables
        while has(equation_final{i}, species)
    
            % Get all the variables in the equation
            all_var = symvar(equation_final{i});
    
            % Get only the steady state species
            dep_var = intersect(all_var, species);
    
            % Replace each steady state species with its corresponding expression
            for j = 1:length(dep_var)
                equation_final{i} = simplify(subs(equation_final{i}, {dep_var(j)}, equation_final{find(species == dep_var(j))}));

                % Check if the dependent variable also has another dependent variable
                dep_var_eq = equation_final{find(species == dep_var(j))};
                all_var2 = symvar(dep_var_eq);
                dep_var2 = intersect(all_var2, species);

                % If it has further dependent variables
                if ~isempty(dep_var2)

                    % Check if the variable to be added is already in the list
                    if ~ismember(dep_var(j), dep_var_cycle)
                        dep_var_cycle{end+1} = dep_var(j);
                    else
                        fprintf('Could not write parametrization in terms of free parameters.\n')
                        fprintf('At least %s and %s are dependent on each other.\n\n', dep_var(j), dep_var_cycle{end})

                        % Indicate that displaying the parametrization needs to be skipped
                        skip_simplification = 1;
                    end
                end
            end

            % Do not continue the while loop if flag is 1
            if skip_simplification == 1

                % Exit the loop
                break
            end
        end

        % Do not continue to the next equation if flag is 1
        if skip_simplification == 1
            break
        end
    end

    % Sort the species
    [species, index] = sort(species);
        
    % Get the free parameters
    free_parameter = string(setdiff(cell2sym(model.species), species));
    
    % If simplification was skipped
    if skip_simplification == 1
    
        % Sort the untouched equations according to the species
        equation = equation(index);

        % Display the untouched steady state solution of the system
        fprintf('The solution is:\n\n')
        for i = 1:length(equation)
            fprintf('%s = %s\n', species(i), equation{i})
        end
    else

        % Otherwise, use the parametrization in terms of free parameters
        equation_final = equation_final(index);

        % Display the parametrized steady state solution of the system
        fprintf('The solution is:\n\n')
        for i = 1:length(equation_final)
            fprintf('%s = %s\n', species(i), equation_final{i})
        end
    end

    % Display the list of free parameters
    if isempty(free_parameter)
        fprintf('Free parameters: None \n\n')
    elseif length(free_parameter) == 1
        fprintf(1, 'Free parameter: ')
        for i = 1:length(free_parameter)-1
            fprintf(1, '%s, ', char(free_parameter(i)'))
        end
        fprintf(1, '%s\n\n', char(free_parameter(end)))
    else
        fprintf(1, 'Free parameters: ')
        for i = 1:length(free_parameter)-1
            fprintf(1, '%s, ', char(free_parameter(i)'))
        end
        fprintf(1, '%s\n\n', char(free_parameter(end)))
    end

    % Display conservation laws
    if isempty(conservation_law_matrix)
        fprintf('Conservation laws: None \n\n')
    elseif length(conservation_law) == 1
        fprintf('Conservation law: %s \n\n', conservation_law{1})
    else
        fprintf('Conservation laws: \n')
        for i = 1:length(conservation_law)
            disp(conservation_law{i})
        end
    end

end










% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                                                   % %
% % The following are functions used in the algorithm % %
% %                                                   % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Function 1 of 13: indepDecomp                                           %
%                                                                         %
%    - Purpose: To decompose a network into its finest independent        %
%         decomposition                                                   %
%    - Input: model: empty species list                                   %
%    - Outputs                                                            %
%         - model: completed structure                                    %
%         - R: matrix of reaction vectors of the network                  %
%         - G: undirected graph of R                                      %
%         - P: partitions representing the decomposition of the reactions %
%    - Used in steadyState (Step 1)                                       %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [model, N, R, G, P] = indepDecomp(model)

    % Create a list of all species indicated in the reactions
    [model, m] = modelSpecies(model);
    
    % Form stoichiometric matrix N
    [N, ~, ~, r] = stoichMatrix(model, m);

    % Get the transpose of N: Each row now represents the reaction vector a reaction    
    R = N';
    
    % Write R in reduced row echelon form: the transpose of R is used so 'basis_reaction_num' will give the pivot rows of R
    %    - 'basis_reaction_num' gives the row numbers of R which form a basis for the rowspace of R
    [~, basis_reaction_num] = rref(R');
    
    % Form a basis for the rowspace of R
    basis = R(basis_reaction_num, :);
    
    % Initialize an undirected graph G
    G = graph();
    
    % Construct the vertex set of undirected graph G
    % Add vertices to G: these are the reaction vectors that form a basis for the rowspace of R
    for i = 1:numel(basis_reaction_num)
        
        % Use the reaction number as label for each vertex
        G = addnode(G, strcat('R', num2str(basis_reaction_num(i))));
    end
    
    % Initialize matrix of linear combinations
    linear_combo = zeros(r, numel(basis_reaction_num));
    
    % Write the nonbasis reaction vectors as a linear combination of the basis vectors
    % Do this for the nonbasis reactions vectors
    for i = 1:r
        if ~ismember(i, basis_reaction_num)
          
          % This gives the coefficients of the linear combinations
          % The basis vectors will have a row of zeros
          linear_combo(i, :) = basis'\R(i, :)';
        end
    end
    
    % Round off values to nearest whole number to avoid round off errors
    linear_combo = round(linear_combo);
        
    % Get the reactions that are linear combinations of at least 2 basis reactions
    % These are the reactions where we'll get the edges
    get_edges = find(sum(abs(linear_combo), 2) > 1);
        
    % Initialize an array for sets of vertices that will form the edges
    vertex_set = { };
     
    % Identify which vertices form edges in each reaction: get those with non-zero coefficients in the linear combinations
    for i = 1:numel(get_edges)
        vertex_set{i} = find(linear_combo(get_edges(i), :) ~= 0);
    end
    
    % Initialize the edge set
    edges = [ ];
    
    % Get all possible combinations (not permutations) of the reactions involved in the linear combinations
    for i = 1:numel(vertex_set)
        edges = [edges; nchoosek(vertex_set{i}, 2)];
    end
    
    % Get just the unique edges
    edges = unique(edges, 'rows');
    
    % Add these edges to graph G
    for i = 1:size(edges, 1)
        G = addedge(G, strcat('R', num2str(basis_reaction_num(edges(i, 1)))), strcat('R', num2str(basis_reaction_num(edges(i, 2)))));
    end
        
    % Determine to which component each vertex belongs to
    component_numbers = conncomp(G);
    
    % Determine the number of connected components of G: this is the number of partitions R will be decomposed to
    num_components = max(component_numbers);
    
    % For the case of only one connected component
    if num_components == 1
        P = [ ];
    end
    
    % Initialize the list of partitions
    P = cell(1, num_components);
    
    % Basis vectors: assign them first into their respective partition based on their component number
    for i = 1:numel(component_numbers)
        P{component_numbers(i)}(end+1) = basis_reaction_num(i);
    end
    
    % Nonbasis vectors: they go to the same partition as the basis vectors that form their linear combination
    for i = 1:numel(P)
        for j = 1:numel(P{i})
            
            % Get the column number representing the basis vectors in 'linear_combo'
            col = find(basis_reaction_num == P{i}(j));
            
            % Check which reactions used a particular basis vector and assign them to their respective partition
            P{i} = [P{i} find(linear_combo(:, col) ~= 0)'];
        end
    end
    
    % Get only unique elements in each partition
    for i = 1:numel(P)
        P{i} = unique(P{i});
    end

    % If some reactions are missing, then redo starting from the linear combination part
    if length(cell2mat(P)) ~= size(R, 1)

        % Do not round off the coefficients of the linear combinations
        linear_combination = linear_combo;

        get_edges = find(sum(abs(linear_combination), 2) > 1);
        vertex_set = { };
        for i = 1:numel(get_edges)
            vertex_set{i} = find(linear_combination(get_edges(i), :) ~= 0);
        end
        edges = [ ];
        for i = 1:numel(vertex_set)
            edges = [edges; nchoosek(vertex_set{i}, 2)];
        end
        edges = unique(edges, 'rows');
        for i = 1:size(edges, 1)
            G = addedge(G, strcat('R', num2str(basis_reaction_num(edges(i, 1)))), strcat('R', num2str(basis_reaction_num(edges(i, 2)))));
        end
        
        component_numbers = conncomp(G);
        num_components = max(component_numbers);
        if num_components == 1
            P = [ ];
        end
        
        P = cell(1, num_components);
        for i = 1:numel(component_numbers)
            P{component_numbers(i)}(end+1) = basis_reaction_num(i);
        end
        for i = 1:numel(P)
            for j = 1:numel(P{i})
                col = find(basis_reaction_num == P{i}(j));
                P{i} = [P{i} find(linear_combination(:, col) ~= 0)'];
            end
        end
        for i = 1:numel(P)
            P{i} = unique(P{i});
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 2 of 13: analyticSolution                                        %
%                                                                           %
%    - Purpose: To solve for the parametrized steady state solution of a    %
%         network                                                           %
%    - Inputs                                                               %
%         - model: empty species list                                       %
%         - P_: list of reaction numbers from the independent decomposition %
%         - B_: list of tau's already used in parametrization               %
%         - sigma_: list of sigma's already used in parametrization         %
%    - Outputs                                                              %
%         - param: list of parametrization of the species of the system     %
%         - B_: updated list of tau's already used in parametrization       %
%         - sigma_: updated list of sigma's already used in parametrization %
%         - k: list of rate constants of the system                         %
%         - sigma: list of sigma's of the system                            %
%         - tau: list of tau's of the system                                %
%         - model: completed structure                                      %
%         - skip: logical; indicator for steadyState to skip solving the    %
%              subnetwork                                                   %
%    - Used in steadyState (Step 2)                                         %
%    - Note: The function uses the class graph_.m (which uses edge.m and    %
%         vertex.m)                                                         %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [param, B_, sigma_, k, sigma, tau, model, skip] = analyticSolution(model, P_, B_, sigma_)

    % Create a list of all species indicated in the reactions
    [model, m] = modelSpecies(model);
    
    % Form stoichiometric matrix N
    [N, reactant_complex, product_complex, r] = stoichMatrix(model, m);

    % Get just the unique complexes
    % index(i) is the index in all_complex of the reactant complex in reaction i
    all_complex = unique([reactant_complex product_complex]', 'rows');
    
    % Construct the matrix of complexes
    all_complex = all_complex';
    
    % Count the number of complexes
    n = size(all_complex, 2);

    % Determine the number of linkage and strong linkage classes
    [l, sl] = linkageClass(reactant_complex, product_complex);
    
    % Check if weakly reversible
    is_weakly_reversible = sl == l;
    
    % Get the rank of the reaction network
    s = rank(N);
    
    % Compute the deficiency of the reaction network
    deficiency = n - l - s;
    
    % Initialize indicator if steadyStage needs to skip the subnetwork
    skip = 0;

    % If the network has deficiency 0 and is weakly reversible, no need to translate
    if is_weakly_reversible & deficiency == 0
        
        % Get the matrix of reactant and product complexes
        stoichiometric_complex_reactant = reactant_complex;
        stoichiometric_complex_product = product_complex;
    
        % Create generalized chemical reaction network (GCRN) vertices
        GCRN_reactant = [stoichiometric_complex_reactant; stoichiometric_complex_reactant];
        GCRN_product = [stoichiometric_complex_product; stoichiometric_complex_product];
        GCRN_vertex = unique([GCRN_reactant, GCRN_product]', 'rows');
        GCRN_vertex = GCRN_vertex';
    else
    
        % Generate the GCRN
        [stoichiometric_complex_reactant, stoichiometric_complex_product, kinetic_complex_reactant, kinetic_complex_product, model, skip] = GCRN(model);

        % If the subnetwork needs to be skipped
        if skip == 1
            
            % Still form rate constants k
            % Create rate constants based on reaction number of decomposition
            for i = 1:length(P_)
                try
                    k(end+1) = sym(strcat('k', string(P_(i))));
                catch
                    k = sym(strcat('k', string(P_(i))));
                end
            end
            
            % Empty results since no translation was found
            sigma = [ ];
            tau = [ ];
            param = [ ];

            % Exit the function
            return
        end

        % Create GCRN vertices
        GCRN_reactant = [stoichiometric_complex_reactant; kinetic_complex_reactant];
        GCRN_product = [stoichiometric_complex_product; kinetic_complex_product];
        GCRN_vertex = unique([GCRN_reactant, GCRN_product]', 'rows');
        GCRN_vertex = GCRN_vertex';
    end
    
    % Initialize graph edges
    GCRN_graph = zeros(2, size(GCRN_reactant, 2));
    
    for i = 1:size(GCRN_reactant, 2)
    
        % Fill out the reactants
        GCRN_graph(1, i) = find(ismember(GCRN_vertex', GCRN_reactant(:, i)', 'rows'));
    
        % Fill out the products
        GCRN_graph(2, i) = find(ismember(GCRN_vertex', GCRN_product(:, i)', 'rows'));
    end
    
    % Initialize spanning forest
    M = zeros(size(stoichiometric_complex_reactant, 1), size(GCRN_vertex, 2) - 1);
    
    % Create spanning forest using kinetic complexes
    for i = 2:size(GCRN_vertex, 2)
       M(:, i-1) = GCRN_vertex(size(GCRN_vertex, 1)/2+1:size(GCRN_vertex, 1), i) - GCRN_vertex(size(GCRN_vertex, 1)/2+1:size(GCRN_vertex, 1), 1);
    end
    
    % Let Mp be M' 
    Mp = M';
    
    % Create a symbolic matrix H
    syms H [size(M, 1) size(M, 2)] matrix
    H = symmatrix2sym(H); % To see the elements of H
    
    % We want to solve H in M' * H * M' = M'
    % LHS
    H_ = M' * H * M';
    
    % Solve for H when LHS = RHS
    H_ = solve([H_(:) == Mp(:)]);
    
    % Extract the solution and place them in matrix H
    % H is the generalized inverse of M'
    for i = 1:length(fieldnames(H_))
        H(i) = getfield(H_, string(H(i)));
    end
    
    % For the case when H has only 1 entry
    if isempty(H(i))
        H = 1;
    end
    
    % Form B: outputs the answer as rational numbers
    % im B = ker M'
    B = null(M', 'r');

    % Number of vertices
    V = size(GCRN_vertex, 2);
    
    % Initialize digraph
    G = graph_(V);
    
    % Initialize edges of spanning trees found
    % The index correspond to the root vertex
    spanTree_edge = cell(1, V);
    
    % Add digraph edges
    for i = 1:size(GCRN_graph, 2)
        G.addEdge(GCRN_graph(1, i), GCRN_graph(2, i));
    end
    
    % Initialize list of tree constant reactions
    tree_constant_reaction = cell(1, V);

    % Get edges of the spanning trees towards each vertex
    for i = 1:V
    
        % Get the spanning trees TOWARDS each vertex
        spanTree = directedSpanTreeTowards(G, i);
    
        % Collect the spanning trees
        for j = 1:length(spanTree)
            spanTree_edge{i}{end+1} = spanTree{j};
        end
        
        % Get reaction numbers for the tree constants
        % Find the corresponding reaction number of the edges
        for j = 1:length(spanTree_edge{i})
            tree_constant_reaction{i}{end+1} = find(ismember(GCRN_graph', spanTree_edge{i}{j}', 'rows'));
        end
    end

    % k: rate constants
    % Create rate constants based on reaction number of decomposition
    for i = 1:length(P_)
        try
            k(end+1) = sym(strcat('k', string(P_(i))));
        catch
            k = sym(strcat('k', string(P_(i))));
        end
    end
    
    % sigma: phantom edge rate constants
    if size(GCRN_graph, 2) - r == 0
        sigma = [ ];
    else
        for i = 1:size(GCRN_graph, 2) - r
        
            % Create tracker of sizes of sigma
            sigma_(end+1) = length(sigma_) + 1;
        
            % Create the sigma's, continuing the numbering from previous subnetwork
            try
                sigma(end+1) = sym(strcat('sigma', string(length(sigma_))));
            catch
                sigma = sym(strcat('sigma', string(length(sigma_))));
            end
        end
    end
    
    % K: tree constants
    K = sym('K', [1 V]);
    
    % kappa
    kappa = sym('kappa', [1 V-1]);
    
    % tau
    if isempty(B)
        tau = [ ];
    else
        for i = 1:size(B, 2)
        
            % Create tracker of sizes of B
            B_(end+1) = length(B_) + 1;
        
            % Create the tau's, continuing the numbering from previous subnetwork
            try
                tau(end+1) = sym(strcat('tau', string(length(B_))));
            catch
                tau = sym(strcat('tau', string(length(B_))));
            end
        end
    end
    
    % Create list of rate constants including for phantom edges (if they exist)
    rate_constants = [k, sigma];
    
    % Form the tree constants
    for i = 1:length(K)
        for j = 1:length(tree_constant_reaction{i})
            if j == 1
                K(i) = prod(rate_constants(tree_constant_reaction{i}{j}));
            else
                K(i) = K(i) + prod(rate_constants(tree_constant_reaction{i}{j}));
            end
            
        end
    end
    
    % Create kappa
    j = 2;
    for i = 1:length(kappa)
        kappa(i) = simplify(K(j)/K(1));
        j = j + 1;
    end
    
    % Generate the parametrization of the steady state
    try
        param = [(kappa.^H.').' (tau.^B.').'];
        param = prod(param.');
    
    % If there are no tau's
    catch
        param = [(kappa.^H.').'];
        param = prod(param.');
    end
    
    % Display the result
    for i = 1:length(model.species)
        fprintf('%s = %s \n', model.species{i}, string(param(i)))
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 3 of 13: modelSpecies                                      %
%                                                                     %
%    - Purpose: To fill out list of species based on given reactions  %
%    - Input: model: empty species list                               %
%    - Outputs                                                        %
%         - model: completed structure                                %
%         - m: number of species                                      %
%    - Used in                                                        %
%         - indepDecomp                                               %
%         - analyticSolution                                          %
%         - GCRN                                                      %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [model, m] = modelSpecies(model)

    % Initialize list of species
    model.species = { };

    % Get all species from reactants
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).reactant)
            model.species{end+1} = model.reaction(i).reactant(j).species;
        end
    end
    
    % Get species from products
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).product)
            model.species{end+1} = model.reaction(i).product(j).species;
        end
    end
    
    % Get only unique species
    model.species = unique(model.species);
    
    % Count the number of species
    m = numel(model.species);
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 4 of 13: stoichMatrix                                            %
%                                                                           %
%    - Purpose: To form the stoichometrix matrix N and the set of reactant  %
%          and product complexes                                            %
%    - Inputs                                                               %
%         - model: complete structure                                       %
%         - m: number of species                                            %
%    - Outputs                                                              %
%         - N: stoichiometric matrix                                        %
%         - reactant_complex: matrix of reactant complexes                  %
%         - product_complex: matrix of product complexes                    %
%         - r: total number of reactions                                    %
%    - Used in                                                              %
%         - indepDecomp                                                     %
%         - analyticSolution                                                %
%         - GCRN                                                            %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [N, reactant_complex, product_complex, r] = stoichMatrix(model, m)

    % Initialize the matrix of reactant complexes
    reactant_complex = [ ];
    
    % Initialize the matrix of product complexes
    product_complex = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
      
        % Initialize the vector for the reaction's reactant complex
        reactant_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complex(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complex(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complex(:, end) - reactant_complex(:, end);
        
        % If the reaction is reversible
        if model.reaction(i).reversible
          
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complex(:, end+1) = product_complex(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complex(:, end+1) = reactant_complex(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Count the total number of reactions
    r = size(N, 2);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 5 of 13: GCRN                                                    %
%                                                                           %
%    - Purpose: To create the GCRN of a translated network                  %
%    - Input: model: empty species list                                     %
%    - Outputs                                                              %
%         - stoichiometric_complex_reactant: reactant stoichiometric        %
%              complex of the GCRN                                          %
%         - stoichiometric_complex_product: product stoichiometric complex  %
%              of the GCRN                                                  %
%         - kinetic_complex_reactant: kinetic complex associated with the   %
%              reactant stoichiometric complex                              %
%         - kinetic_complex_product: kinetic complex associated with the    %
%              product stoichiometric complex                               %
%         - model: completed structure                                      %
%         - skip: logical; indicator for steadyState to skip solving the    %
%              subnetwork                                                   %
%    - Used in analyticSolution                                             %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [stoichiometric_complex_reactant, stoichiometric_complex_product, kinetic_complex_reactant, kinetic_complex_product, model, skip] = GCRN(model)

    % Create a list of all species indicated in the reactions
    [model, m] = modelSpecies(model);
    
    % Form stoichiometric matrix N
    [~, reactant_complex, product_complex, r] = stoichMatrix(model, m);

    % Get the translations of the network
    [Solution, Index, skip] = findTranslations(reactant_complex, product_complex);

    % If the subnetwork needs to be skipped
    if skip == 1

        % Empty results since no translation was found
        stoichiometric_complex_reactant = [ ];
        stoichiometric_complex_product = [ ];
        kinetic_complex_reactant = [ ];
        kinetic_complex_product = [ ];
        return
    end

    % Go through each translation
    for i = 1:size(Solution, 1)

        % Get the translation
        translated_network = Solution(i, :);

        % Get the Index
        Index_ = Index(i, :);

        % Split merged reactions
        % Initialize expanded reactant and product matrices of the translated network
        translated_network_reactant_expanded = [ ];
        translated_network_product_expanded = [ ];

        % Split merged reactions
        for j = 1:length(Index_)
            if length(Index_{j}) == 1
                translated_network_reactant_expanded(:, end+1) = translated_network{1}(:, j);
                translated_network_product_expanded(:, end+1) = translated_network{2}(:, j);
            elseif length(Index_{j}) > 1
                for k = 1:length(Index_{j})
                    translated_network_reactant_expanded(:, end+1) = translated_network{1}(:, j);
                    translated_network_product_expanded(:, end+1) = translated_network{2}(:, j);
                end
            end
        end
    
        % Get the index which will be used to rearrange the complexes
        index_ = cell2mat(Index_);
        
        % Rearrange reactant and product complexes of the translation
        % Place the complexes in correct order of reactions
        [~, rearrange] = ismember(1:numel(index_), index_);
        translated_network_reactant = translated_network_reactant_expanded(:, rearrange);
        translated_network_product = translated_network_product_expanded(:, rearrange);

        % Get unique complexes of the translation
        translation_complex = unique([translated_network_reactant, translated_network_product]', 'rows');
        translation_complex = translation_complex';
    
        % Initialize list of count of kinetic complexes per stoichiometric complex
        % The index corresponds to the complex in translation_complex
        kinetic_count = zeros(1, size(translation_complex, 2));
    
        % Initialize list of translation reactant complexes which occurs multiple times
        translation_reactant_complex_count = zeros(1, size(translation_complex, 2));
        
        % Assign kinetic complexes to each stoichiometric complex
        % Go through each unique stoichiometric complex
        for j = 1:numel(kinetic_count)
    
            % Look for the index of this complex in translated_network_reactant
            index_in_translation_reactant = find(ismember(translated_network_reactant', translation_complex(:, j)', 'rows'));
    
            % Assign the number of indices to translation_reactant_complex_multiple
            translation_reactant_complex_count(j) = numel(index_in_translation_reactant);
    
            % Assign the kinetic complex of these reaction numbers to the stoichiometric complex
            kinetic = unique(reactant_complex(:, index_in_translation_reactant)', 'rows')';
    
            % Count how many kinetic complexes
            kinetic_count(j) = size(kinetic, 2);
        end
    
        % Construct a list of reactant complexes in the translation that appear multiple times
        translation_reactant_complex_multiple = find(translation_reactant_complex_count > 1);
    
        % Check if a stoichiometric complex has more than one kinetic complex
        % If yes, no need to add phantom edges
        if ~any(kinetic_count > 1)

            % Form matrices of stoichiometric complexes
            stoichiometric_complex_reactant = translated_network_reactant;
            stoichiometric_complex_product = translated_network_product;
    
            % Form matrices of kinetic complexes
            kinetic_complex_reactant = reactant_complex;
            kinetic_complex_product = zeros(m, r);
    
            % For the matrix of product kinetic complexes, create a matrix of vertex numbers matching reactants and products
            [~, reactant_vertex] = ismember(translated_network_reactant', translation_complex', 'rows');
            [~, product_vertex] = ismember(translated_network_product', translation_complex', 'rows');
            translation_vertex = [reactant_vertex'; product_vertex'];
    
            % Initialize GCRN vertices
            GCRN_vertices = zeros(2, r);
    
            % Fill out the vertex numbers for the reactant
            GCRN_vertices(1, :) = 1:r;
    
            % Fill out the vertex numbers for the product: match them with the numbering of the row for reactants
            for j = 1:r
    
                % Simply get the first in the list of reaction numbers if there are multiple ones
                if ismember(translation_vertex(2, j), translation_reactant_complex_multiple)
                    reaction_number = find(translation_vertex(1, :) == translation_vertex(2, j));
                    GCRN_vertices(2, j) = reaction_number(1);
                else
                    GCRN_vertices(2, j) = find(translation_vertex(1, :) == translation_vertex(2, j));
                end
            end
            
            % Finish the matrix of product kinetic complexes
            kinetic_complex_product(:, 1:r) = kinetic_complex_reactant(:, GCRN_vertices(2, 1:r));
    
            % Compute the kinetic deficiency
            kinetic_deficiency = kineticDeficiency(kinetic_complex_reactant, kinetic_complex_product);
    
            % End once the kinetic deficiency is 0
            if kinetic_deficiency == 0
                break
            end
    
        % Otherwise, phantom edges are needed
        else

            % Form matrices of stoichiometric complexes
            stoichiometric_complex_reactant = translated_network_reactant;
            stoichiometric_complex_product = translated_network_product;
    
            % Form matrices of kinetic complexes
            kinetic_complex_reactant = reactant_complex;
            kinetic_complex_product = zeros(m, r);
    
            % For the matrix of product kinetic complexes, create a matrix of vertex numbers matching reactants and products
            [~, reactant_vertex] = ismember(translated_network_reactant', translation_complex', 'rows');
            [~, product_vertex] = ismember(translated_network_product', translation_complex', 'rows');
            translation_vertex = [reactant_vertex'; product_vertex'];
    
            % Initialize GCRN vertices
            GCRN_vertices = zeros(2, r);
    
            % Fill out the vertex numbers for the reactant
            GCRN_vertices(1, :) = 1:r;
    
            % Determine the complexes that need phantom edges
            complex_need_phantom_edge = find(kinetic_count > 1);
    
            % Fill out the vertex numbers for the product: match them with the numbering of the row for reactants
            for j = 1:r
    
                % If the complex needs a phantom edge, initialize it first to 0
                if ismember(translation_vertex(2, j), complex_need_phantom_edge)
                    GCRN_vertices(2, j) = 0;
                else

                    % This ensures that we get the first in the list of reaction numbers if there are multiple ones but doesn't need a phantom edge
                    reaction_number = find(translation_vertex(1, :) == translation_vertex(2, j));
                    GCRN_vertices(2, j) = reaction_number(1);
                end
            end
            
            % Go through each complex that needs phantom edges
            for j = 1:numel(complex_need_phantom_edge)
            
                % Get the complex's list of corresponding reaction numbers
                list_reaction_number = find(translation_vertex(1, :) == complex_need_phantom_edge(j));
            
                % Get first reaction: this will be kept
                keep_index = list_reaction_number(1);
            
                % Remove this reaction from the list so that the remaining ones are those that will need phantom edges
                list_reaction_number(1) = [ ];
            
                % Add phantom edges
                % Go through the rest of the reactions left
                for k = 1:numel(list_reaction_number)
            
                    % Get the reaction number
                    % This reaction is deleted, then a phantom edge is constructed, then the reaction is returned
                    delete_index = list_reaction_number(k);
            
                    % Create phantom edge
                    GCRN_vertices(1, end+1) = size(GCRN_vertices, 2) + 1;
                    GCRN_vertices(2, end) = delete_index;
            
                    % Updated the matrices of stoichiometric and kinetic complexes
                    stoichiometric_complex_reactant(:, end+1) = stoichiometric_complex_reactant(:, keep_index);
                    stoichiometric_complex_product(:, end+1) = stoichiometric_complex_reactant(:, keep_index);
                    kinetic_complex_reactant(:, end+1) = kinetic_complex_reactant(:, keep_index);
                    kinetic_complex_product(:, end+1) = kinetic_complex_reactant(:, delete_index);
                end
            
                % Look for the complex in the second row of translation_vertex
                need_replacement = find(translation_vertex(2, :) == complex_need_phantom_edge(j));
            
                % For those we kept at 0 earlier
                for k = 1:numel(need_replacement)
    
                    % If the corresponding vertex was not used to create a phantom edge, we use the reaction that was kept
                    if numel(find(GCRN_vertices(1, :) == need_replacement(k))) == 1
                        GCRN_vertices(2, need_replacement(k)) = keep_index;
    
                    % Otherwise, we use the reaction removed then replaced back
                    else
                        GCRN_vertices(2, need_replacement(k)) = delete_index;
                    end
                end
            end
            
            % Finish the matrix of product kinetic complexes
            kinetic_complex_product(:, 1:r) = kinetic_complex_reactant(:, GCRN_vertices(2, 1:r));
    
            % Compute the kinetic deficiency
            kinetic_deficiency = kineticDeficiency(kinetic_complex_reactant, kinetic_complex_product);
    
            % End once the kinetic deficiency is 0
            if kinetic_deficiency == 0
                break
            end
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 6 of 13: findTranslations                                        %
%                                                                           %
%    - Purpose: To determine the different weakly reversible and deficiency %
%         zero translations of a CRN                                        %
%    - Inputs                                                               %
%         - reactant_complex: matrix of reactant complexes                  %
%         - product_complex: matrix of product complexes                    %
%    - Outputs                                                              %
%         - Solution: reactant complexes and product complexes of each      %
%              translation                                                  %
%         - Index: permutation of the reaction numbers between the original %
%              CRN and the translation                                      %
%         - skip: logical; indicator for steadyState to skip solving the    %
%              subnetwork                                                   %
%    - Used in GCRN                                                         %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [Solution, Index, skip] = findTranslations(reactant_complex, product_complex)

    % Count the number of species
    m = size(reactant_complex, 1);
    
    % Form the stroichiometric matrix
    N = product_complex - reactant_complex;
    
    % Determine the rank of the CRN
    s = rank(N);
    
    % Initialize order of complexes
    max_order = 2;

    % Keep on increasing max_order until 'Solution' is non-empty
    while true

        % Prepare a matrix of all possible combinations of m elements out of a list from 1 to m+max_order
        tmp_mat =  nchoosek(1:(m+max_order), m)';
        
        % Prepare a matrix of NaN's of size m x nchoosek(m+max_order, max_order)
        total_complexes = nan(m, nchoosek(m+max_order, max_order));
        
        % Get the number of complexes in total_complexes
        cmplx_num = size(total_complexes, 2);
        
        % Generate all posible complexes under a given maximum order, 'max_order'
        for j = 1:cmplx_num
            total_complexes(1, j) =  tmp_mat(1, j) - 1;
            for i = 2:m
                total_complexes(i, j) =  tmp_mat(i, j) - tmp_mat(i-1, j) - 1;
            end
        end
        
        % Prepare a matrix of NaN's of size m x cmplx_num*(cmplx_num - 1)/2
        tmp_mat = nan(m, cmplx_num*(cmplx_num - 1)/2);
        
        % Create all the stoichiometric vectors for all possible complexes
        cid = 0;
        for i = 1:cmplx_num
            for j = (i+1):cmplx_num
                cid = cid+1;
                tmp_mat(:, cid) = total_complexes(:, j) - total_complexes(:, i);
            end
        end
        
        % 'base_edge' contains all the stoichiometric vectors for all possible complexes; however, it contains only either -v or v, but not both
        base_edge = unique(tmp_mat', 'rows')';
        
        % Random number generator at seed 2
        rng(2);
    
        % If the network is NOT empty
        if ~isempty(reactant_complex) == 1
        
            % Extend N to include its reverse
            N_extend = [N, -N];
        
            % Check the necessary condition for having zero deficiency after translation (Theorem 3.3 and 3.5)
            base_edge_id = ismember(base_edge', N_extend', 'rows');
            comp_reject  = 0;
            if sum(base_edge_id) > m * (2*m - 1) % Theorem 3.5
                comp_reject = 1;
            elseif sum(base_edge_id) > 2 * m % Theorem 3.3
                comp_reject = 1;
                for ss = 3:(2*m)
                    test_list = nchoosek(find(base_edge_id), ss);
                    sign_matrix = ones(2^(ss-1), ss-1);
                    for rr = 0:(2^(ss-1) - 1)
                        sign_matrix(rr+1, :) = 1 -  2 * decimalToBinary(rr, ss-1);
                    end
                    
                    for rr = 1:size(test_list, 1)
                        base_edge_test_id = test_list(rr, :);
                        base_edge_test = base_edge(:, base_edge_test_id);
                        for sg = 1:2^(ss-1)
                            base_edge_test_sign = repmat([1, sign_matrix(sg, :)], [m, 1]) .* base_edge_test;
                            if isequal(sum(base_edge_test_sign, 2), zeros(m, 1))
                                comp_reject = 0;
                                break
                            end
                        end
                        if comp_reject == 0
                            break
                        end
                    end
                    if comp_reject == 0
                        break
                    end
                end
            end
            
            % Initialize Solution: Contains the source and product complex matrices of found translated CRNs with WR and ZD
            Solution = { };
        
            % Initialize Index: Contains the permutation of reactions between the original CRN and the found CRNs
            Index = { };
        
            if size(N, 2) == s & comp_reject
                Solution = { };
                Index = { };
            else
                
                % Checks the necessary condition for having weak reversibility after translation (Theorem 3.1)
                [~, pos_sol_TF] = positiveSolution(N);
                
                if pos_sol_TF == 0 & comp_reject
                    Solution = { };
                    Index = { };
                else
                    
                    % Perform network translation (Step 2 & 3)
                    [Solution, Index, skip] = translation(reactant_complex, product_complex, max_order);

                    % If the subnetwork needs to be skipped
                    if skip == 1

                        % Exit the function
                        return
                    end
                end
            end
        end

        % Stop looking when a solution is already found
        if ~isempty(Solution)
            break

        % Otherwise, increase the maximum order
        else
            max_order = max_order + 1;
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 7 of 13: kineticDeficiency                                       %
%                                                                           %
%    - Purpose: To determine the kinetic deficiency of a GCRN               %
%    - Inputs                                                               %
%         - kinetic_complex_reactant: matrix of kinetic reactant complexes  %
%         - kinetic_complex_product: matrix of kinetic product complexes    %
%    - Output: kinetic_deficiency: kinetic deficiency                       %
%    - Used in GCRN                                                         %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function kinetic_deficiency = kineticDeficiency(kinetic_complex_reactant, kinetic_complex_product)

    % Get unique kinetic complexes
    complex_translated = unique([kinetic_complex_reactant, kinetic_complex_product]', 'rows');
    complex_translated = complex_translated';
    
    % Get number of kinetic complexes
    n_translated = size(complex_translated, 2);
    
    % Get number of linkage classes
    [l_translated, ~] = linkageClass(kinetic_complex_reactant, kinetic_complex_product);
    
    % Get rank of the network of kinetic complexes
    s_translated = rank(kinetic_complex_product - kinetic_complex_reactant);
    
    % Compute the kinetic deficiency
    kinetic_deficiency = n_translated - l_translated - s_translated;

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 8 of 13: linkageClass                                      %
%                                                                     %
%    - Purpose: To determine the number of linkage and strong linkage %
%         classes                                                     %
%    - Inputs                                                         %
%         - reactant_complex: matrix of reactant complexes            %
%         - product_complex: matrix of product complexes              %
%    - Outputs                                                        %
%         - l: number of linkage classes                              %
%         - sl: number of strong linkage classes                      %
%    - Used in                                                        %
%         - analyticSolution                                          %
%         - kineticDeficiency                                         %
%         - translation                                               %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [l, sl] = linkageClass(reactant_complex, product_complex)

    % Get the number of reactions
    r = size(reactant_complex, 2); 
    
    % Get just the unique complexes
    % index(i) is the index in 'complex' of the reactant complex in reaction i
    [complex, ~, index] = unique([reactant_complex, product_complex]', 'rows');
    
    % Construct the matrix of complexes
    complex = complex';
    
    % Count the number of complexes
    n = size(complex, 2);

    % Initialize a matrix (complexes x complexes) for the reacts_to relation
    % This is for testing reversibility of the network
    reacts_to = false(n, n);
    
    % Initialize matrix (complexes x total reactions) for the reacts_in relation
    % This is the incidence matrix
    reacts_in = zeros(n, r);
    
    % Fill out the entries of the matrices
    for i = 1:r
        
        % reacts_to(i, j) = true iff there is a reaction r: y_i -> y_j
        reacts_to(index(i), index(i + r)) = true;
        
        % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
        reacts_in(index(i), i) = -1;
        reacts_in(index(i+r), i) = 1;
    end
    
    % Linkage classes
    % Count number of connected components of an undirected graph
    linkage_class = conncomp(graph(reacts_to | reacts_to'));
    
    % Count the number of linkage classes
    l = max(linkage_class);

    % Check if the network is reversible
    is_reversible = isequal(reacts_to, reacts_to');
    
    % Strong linkage classes
    % Count number of strongly connected components of an directed graph
    if is_reversible
        strong_linkage_class = linkage_class;
    else
        % Count number of connected components of a directed graph
        strong_linkage_class = conncomp(digraph(reacts_to));
    end
    
    % Count the number of strong linkage classes
    sl = max(strong_linkage_class);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Function 9 of 13: decimalToBinary                                       %
%                                                                         %
%    - Purpose: To convert a number in the decimal system (i.e., base 10) %
%         to the binary number (i.e., base 2)                             %
%    - Inputs                                                             %
%         - decimal_num: a number in decimal system                       %
%         - bin_len: desired length of binary number                      %
%    - Output: binary_num: the number in binary system                    %
%    - Used in findTranslations                                           %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function binary_num = decimalToBinary(decimal_num, bin_len)

    % Note: 0 <= decimal_num <= 2^bin_len - 1

    % Initialize the binary number
    binary_num = zeros(1, bin_len);    

    num_tmp = decimal_num;
    for i = (bin_len-1):-1:0
        if num_tmp >= 2^i
           binary_num(bin_len - i) = 1;
           num_tmp = num_tmp - 2^i;
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 10 of 13: positiveSolution                                       %
%                                                                           %
%    - Purpose: To check the existence of the vector with positive entries  %
%         in the kernel of a given matrix                                   %
%    - Input: Aeq: a matrix                                                 %
%    - Outputs                                                              %
%         - pos_sol: the positive solution to Aeq * x = 0                   %
%         - pos_sol_TF: logical; whether or not Aeq * x = 0 has a positive  %
%              solution                                                     %
%    - Used in findTranslations                                             %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [pos_sol, pos_sol_TF] = positiveSolution(Aeq)

    % Determine the number of rows and columns of the given matrix
    [numrow, numcol] = size(Aeq);

    % If the stoichiometric matrix has the full-column rank, there is only the trivial solution
    if rank(Aeq) == numcol
        pos_sol_TF = 0;
        pos_sol = 0;
    else
        fun0 = @(x) x;
        x0 = zeros(numcol, 1);
        A = [ ];
        b = [ ];
        lb = -ones(numcol, 1);
        ub = zeros(numcol, 1);
        beq = zeros(numrow, 1);
        nonlcon = [ ];
        options = optimoptions('fminimax', 'Display', 'none');
        x = -fminimax(fun0, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
        
        if isequal(x, zeros(numcol, 1))
            pos_sol_TF = 0;
            pos_sol = 0;
        else
            pos_sol_TF = 1;
            pos_sol = x;
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 11 of 13: translation                                            %
%                                                                           %
%    - Purpose: To generate all possible translated networks and check weak %
%         reversibility and zero deficiency                                 %
%    - Inputs                                                               %
%         - reactant_complex: matrix of reactant complexes                  %
%         - product_complex: matrix of product complexes                    %
%         - max_order: maximum number of species in a complex               %
%    - Outputs                                                              %
%         - Solution: list of reactant complexes and product complexes of   %
%              each translation (represented by each row)                   %
%         - Index: list of permutations of the reaction numbers between the %
%              original CRN and the translation                             %
%         - skip: logical; indicator for steadyState to skip solving the    %
%              subnetwork                                                   %
%    - Used in findTranslations                                             %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [Solution, Index, skip] = translation(reactant_complex, product_complex, max_order)
    
    % m: the number of species
    % r: the number of reactions
    [m, r] = size(reactant_complex);
    
    % Form the stoichiometric matrix
    N = product_complex - reactant_complex;

    % Compute for the rank of the network
    s = rank(N);
    
    % Initialize reduced reactant and product complexes
    reactant_complex_reduced = -ones(m, r);
    product_complex_reduced = -ones(m, r);
    
    % For each reaction, reactant_complex_reduced and product_complex_reduced have disjoint supports because all common species are subtracted
    for i = 1:r
        reactant_complex_reduced(:, i) = reactant_complex(:, i) - min(reactant_complex(:, i), product_complex(:, i));
        product_complex_reduced(:, i) = product_complex(:, i) - min(reactant_complex(:, i), product_complex(:, i));
    end
    
    tmp = nchoosek(1:(m+max_order), m);
    max_copy_number = size(tmp, 1);
    additional_complex_tmp = (diff([zeros(max_copy_number, 1), tmp], 1, 2) - 1)';
    additional_complex = -ones(m, nchoosek(m+max_order-1, m-1), max_order);
    order_number = zeros(max_order, 1);
    for i = 2:size(additional_complex_tmp, 2)
        order_tmp = sum(additional_complex_tmp(:, i));
        order_number(order_tmp) = order_number(order_tmp) + 1;
        additional_complex(:, order_number(order_tmp), order_tmp) = additional_complex_tmp(:, i);
    end
    
    reactant_complex_copy = -ones(m, max_copy_number, r);
    product_complex_copy = -ones(m, max_copy_number, r);
    order_number_cumul = cumsum(order_number);
    reaction_orders = max(sum(reactant_complex_reduced, 1), sum(product_complex_reduced, 1));
    
    % Consider both reactant_complex and product_complex for the reaction orders: it is valid because the resulting network is weakly reversible
    for i = 1:r
        reactant_complex_copy(:, 1, i) = reactant_complex_reduced(:, i);
        product_complex_copy(:, 1, i) = product_complex_reduced(:, i);
        for j = 1:(max_order - reaction_orders(i))
            if j == 1
                reactant_complex_copy(:, 2:(order_number_cumul(j)+1), i) = reactant_complex_reduced(:, i) + additional_complex(:, 1:order_number(j), j);
                product_complex_copy(:, 2:(order_number_cumul(j)+1), i) = product_complex_reduced(:, i) + additional_complex(:, 1:order_number(j), j);
            else
                reactant_complex_copy(:, (order_number_cumul(j-1)+2):(order_number_cumul(j)+1), i) = reactant_complex_reduced(:, i)  + additional_complex(:, 1:order_number(j), j);
                product_complex_copy(:, (order_number_cumul(j-1)+2):(order_number_cumul(j)+1), i) = product_complex_reduced(:, i)  + additional_complex(:, 1:order_number(j), j);
            end
        end
    end
    number_of_copies = zeros(1, r);
    for i = 1:r
        if reaction_orders(i) == max_order
            number_of_copies(i) = 1;
        else
            number_of_copies(i) = order_number_cumul(max_order - reaction_orders(i)) + 1;
        end
    end
    number_of_comb = prod(number_of_copies);
    
    sources_tmp = -ones(m, r);
    products_tmp = -ones(m, r);
    
    sol_idx = 1;
    
    Solution = { };
    Index = { };
    
    current_copy = ones(r, 1);
    current_copy(1) = 0;

    % Time the translation process
    iniTime = clock;

    % Time limit to generate translations (set to 10 seconds)
    limit = 10;

    % Initialize indicator if steadyStage needs to skip the subnetwork
    skip = 0;
    
    % Create a matrix containing a vector between ones(r, 1) and number_of_copies
    % We use the upper bound 10^9 to avoid indefinitely running code
    for i = 1:min(number_of_comb, 10^9) 
        
        current_copy(1) = current_copy(1) + 1;
        for j = 1:r
            if current_copy(j) > number_of_copies(j)
                current_copy(j) = 1;
                current_copy(j+1) = current_copy(j+1) + 1;
            end
        end
        
        for j = 1:r
            sources_tmp(:, j) = reactant_complex_copy(:, current_copy(j), j);
            products_tmp(:, j) = product_complex_copy(:, current_copy(j), j);
        end
        
        n = size(unique([sources_tmp, products_tmp]', 'rows')',2);
        
        % If any of the below conditions holds then the deficiency cannot be 0
        if n < s + 1
            continue
        elseif n > 2*s
            continue
        end
        
        % Determine the number of linkage and strong linkage classes
        [l, sl] = linkageClass(sources_tmp, products_tmp);

        % Compute the deficiency of the network
        deficiency = n - l - s;
        if deficiency == 0 & sl == l
            sol_idx = sol_idx + 1;
            complexes_tmp = [sources_tmp; products_tmp];
            [complexes_unique, ~, ic] = unique(complexes_tmp', 'rows');
            complexes_unique = complexes_unique';
            sources_unique = complexes_unique(1:m, :);
            products_unique = complexes_unique((m+1):(2*m), :);
            K_trans = size(sources_unique, 2);
            Solution = [Solution; {sources_unique, products_unique}];
            indexset = cell(1, K_trans);
            for j = 1:r
                indexset{ic(j)} = [indexset{ic(j)}, j];
            end
            Index = [Index; indexset];
        end

        % Check the time limit before moving to the next i
        if etime(clock, iniTime) < limit

            % Continue to the next iteration if within time limit
            continue
        else

            % Indicate that steadyState needs to skip the subnetwork
            skip = 1;

            % Exit the function
            return
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 12 of 13: directedSpanTreeTowards                                %
%                                                                           %
%    - Purpose: To create a list of directed spanning trees towards the     %
%         indicated root vertex of a directed graph                         %
%    - Inputs                                                               %
%         - G: directed graph with vertices and edges                       %
%         - r: root vertex                                                  %
%    - Output: spanTree: edges of the directed spanning tree towards the    %
%         root vertex                                                       %
%    - Used in analyticSolution                                             %
%    - Note: The function uses the class graph_.m (which uses edge.m and    %
%         vertex.m)                                                         %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function spanTree = directedSpanTreeTowards(G, r)

    % Create a new graph for the reversed edges
    G2 = graph_(G.V);

    % Initialize reversed edges
    reverse_edge = cell(1, G.V);
    
    % Reverse the edges of G
    for i = 1:G.V
        for j = 1:length(G.edge{i})
            reverse_edge{G.edge{i}(j)}(end+1) = i;
        end
    end
    
    % Replace the edges of G
    G2.edge = reverse_edge;

    % Root vertex where the spanning trees will come from
    G2.root_vertex = r;

    % Initialize subtree for spanning tree generation
    T = graph_(G2.V);

    % Initialize latest spanning tree found
    L = graph_(G2.V);

    % Initialize list of all edges directed from vertices in T to vertices not in T
    F = [ ];

    % Initialize count of spanning trees rooted at G.root_vertex
    spanTreeCount = 0;

    % Initialize list of spanning trees found
    spanTree = { };
    
    % Initialize F to contain all edges from the root vertex
    % Go through each vertex to which the root vertex is connected to
    for i = 1:length(G2.edge{G2.root_vertex})
    
        % Insert at the beginning of F each edge from the root vertex
        F = [edge(G2.root_vertex, G2.edge{G2.root_vertex}(i)), F];
    end
    
    % Generate spanning trees
    spanTree = grow(G2.V, G2, T, L, F, spanTreeCount, spanTree);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 13 of 13: grow                                                   %
%                                                                           %
%    - Purpose: To generate the spanning trees away from the root vertex    %
%    - Inputs                                                               %
%         - V: number of vertices in the graph                              %
%         - G: given digraph                                                %
%         - T: subtree for spanning tree generation                         %
%         - L: latest spanning tree found                                   %
%         - F: list of all edges directed from vertices in T to vertices    %
%              not in T                                                     %
%         - spanTreeCount: number of spanning trees                         %
%         - spanTree: list of edges of the directed spanning tree towards   %
%              the root vertex                                              %
%    - Outputs                                                              %
%         - spanTree: list of edges of the directed spanning tree towards   %
%              the root vertex                                              %
%         - spanTreeCount: count of spanning trees rooted at G.root_vertex  %
%    - Used in directedSpanTreeTowards                                      %
%    - Notes                                                                %
%         - The function uses the class edge.m                              %
%         - This function was converted from Python to Matlab               %
%         - Python code can be found in https://github.com/ricelink/        %
%              finding-all-spanning-trees-in-directed-graph                 %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
 
function [spanTree, spanTreeCount] = grow(V, G, T, L, F, spanTreeCount, spanTree)

    % Make sure the root vertex of the subgraphs are the same as that of the main graph
    T.root_vertex = G.root_vertex;
    L.root_vertex = G.root_vertex;

    % Check if subtree T already contains all the vertices
    if T.connectedVertices == V

        % If yes, then T becomes L, the latest spanning tree found
        % We deepcopy T to L
        L.V = T.V;
        L.root_vertex = T.root_vertex;
        L.edge = T.edge;
        L.time = T.time;
        L.vertex = T.vertex;

        % Initialize the reversed back edges
        reverse_edge_back = cell(1, V);
        
        % Reverse back the edges
        for i = 1:V
            for j = 1:length(L.edge{i})
                reverse_edge_back{L.edge{i}(j)}(end+1) = i;
            end
        end
        
        % Replace the edges of L
        L.edge = reverse_edge_back;

        % Add 1 to the count of spanning trees
        spanTreeCount = spanTreeCount + 1;

        % Initialize vertex of spanning tree reactions
        spanTreeReaction = zeros(2, length(cell2mat(L.edge)));
        
        % Get the edges from the spanning tree
        i = 1;
        for j = 1:L.V
            if ~isempty(L.edge{j})
                for k = 1:length(L.edge{j})
                    spanTreeReaction(1, i) = j;
                    spanTreeReaction(2, i) = L.edge{j}(k);
                    i = i + 1;
                end
            end
        end
        
        % Add the spanning tree to the list
        spanTree{spanTreeCount} = spanTreeReaction;
        
        % Initialize the reversed reverse_edge_back
        reverse_edge_back_again = cell(1, V);
        
        % Reverse back the edges
        for i = 1:V
            for j = 1:length(L.edge{i})
                reverse_edge_back_again{L.edge{i}(j)}(end+1) = i;
            end
        end
        
        % Replace the edges of L
        L.edge = reverse_edge_back_again;

        % Visit all the vertices of L
        L.depthFirstSearch(L.root_vertex);
    else

        % Initialize FF, a list of edges
        FF = [ ];

        % Run the following block of codes until it breaks
        while true
                
            % If F is not empty, pop an edge e from F
            if ~isempty(F)
    
                % Get the first element of F (an edge)
                e = F(1);
    
                % Remove the first element from F
                F(1) = [ ];
            else
    
                % If F is empty, then the function ends
                return
            end
    
            % Get the ending node of edge e
            % This vertex v is not yet in T
            v = e.to_node;
    
            % Add edge e to T
            T.addEdge(e.from_node, v);

            % Deepcopy F to F_copy
            F_copy = [ ];
            F_copy = F(:)';
    
            % Push each edge (v, w), where w is NOT in T, onto F
            % Go through each vertex to which vertex v is connected to
            for i = 1:length(G.edge{v})
    
                % Get each ending node
                w = G.edge{v}(i);
    
                % If vertex w is not a vertex of any edge in T
                if T.connectedVertex(w) == 0
    
                    % Insert the edge at the beginning of F
                    F = [edge(v, w), F];
                end
            end
            
            % Remove each edge (w, v), where w IS in T, from F
            % Go through each vertex w
            for w = 1:length(G.edge)
    
                % If vertex w is a vertex of an edge in T
                if T.connectedVertex(w) == 1
    
                    % Initialize checker if w -> v exists in T
                    w_to_v = [ ];
    
                    % Go through each vertex to which vertex w is connected to
                    for j = 1:length(G.edge{w})
    
                        % If vertex w is connected to vertex v
                        if G.edge{w}(j) == v
    
                            % Append (i.e., place at the end) the index of v to the checker
                            w_to_v = [w_to_v, j];
                        end
                    end
    
                    % If there is w -> v in T
                    if ~isempty(w_to_v)
    
                        % Initialize checker of location of w -> v in F
                        w_to_v_loc = [ ];
    
                        % Go through each element of F
                        for j = 1:length(F)
    
                            % Look for the location of w -> v in F
                            if (F(j).from_node == w) & (F(j).to_node == v)
    
                                % Append the index in F to the checker
                                w_to_v_loc = [w_to_v_loc, j];
                            end
                        end
    
                        % If w -> v is located in F
                        if ~isempty(w_to_v_loc)
    
                            % j decreases by 1 each step until it reaches 1
                            for j = length(w_to_v_loc):-1:1
                                
                                % Remove the w -> v from F
                                F(w_to_v_loc(j)) = [ ];
                            end
                        end
                    end
                end
            end
    
            % Recurse
            [spanTree, spanTreeCount] = grow(G.V, G, T, L, F, spanTreeCount, spanTree);

            % Pop each edge (v, w), where w is NOT in T, from F and restore each edge (w, v), where w IS in T, in F
            % Deepcopy F_copy to F
            F = F_copy(:)';

            % Remove edge e from graphs T and G
            T.removeEdge(e.from_node, e.to_node);
            G.removeEdge(e.from_node, e.to_node);

            % Insert edge e at the beginning of FF
            FF = [e, FF];
    
            % Bridge Test
            %    - An edge e is a bridge if G\{e} is NOT connected
            %    - Check if there is an edge (w, v), where w is NOT a descendant of v in L

            % Initialize b as true
            b = 1;
    
            % Go through each vertex
            for w = 1:length(G.edge)
    
                % Variable to check for an edge w -> v
                w_to_v = -1;
    
                % Go through each vertex to which vertex w is connected to
                for j = 1:length(G.edge{w})
    
                    % If vertex w is connected to vertex v
                    if G.edge{w}(j) == v
    
                        % j is the index of v s.t. w -> v
                        % This confirms that w is NOT a descendant of v
                        w_to_v = j;
                    end
                end

                % If w is NOT directly connected to v
                if w_to_v ~= -1

                    % Check if w is a descendant of any child of v in L
                    if ((L.vertex{v}.d < L.vertex{w}.d) & (L.vertex{w}.d < L.vertex{w}.f) & (L.vertex{w}.f < L.vertex{v}.f)) == 0
    
                        % Change b to 0 (i.e., the Bridge Test fails)
                        b = 0;
    
                        % Get out of the for-loop for w
                        break
                    end
                end
            end
    
            % If b is still 1 (i.e., the Bridge Test passes)
            if b == 1

                % Get out of the while-true loop
                break
            end
        end
    
        % Reconstruct G
        while ~isempty(FF)
    
            % Get the first element of FF (an edge)
            e = FF(1);
    
            % Remove it from FF
            FF(1) = [ ];
    
            % Insert it at the beginning of F
            F = [e, F];
    
            % Add it to G
            G.addEdge(e.from_node, e.to_node);
        end
    end

end