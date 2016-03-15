function [Customer2Item,Item2Customer,Nsolutions,Rewards]=Murty(RewardMatrix,N)
%Implementation of Murty's Algorithm to find N-best solutions to an
%assignment problem
%by Umut Orguner

%Requires Auction.m

%Inputs
%-RewardMatrix: Assignment matrix for which the best N solutions are
%required
%-N: The number of best solutions required from the assignment problem
%defined by RewardMatrix

%Outputs
%-Customer2Item: List of solutions similar to Auction.m -- Customer
%Assignments. Each row corresponds to a different solution.
%-Item2Customer: List of solutions similar to Auction.m -- Item Assignments
%(Assignment 0 means no assignment) Each row corresponds to a different
%solution.
%-Nsolutions: The number of solutions (rows) output in Customer2item and Item2Customer.  
%The number of solutions are limited. So if too many solutions are 
%required as input (N), the maximum number of solutions are returned which
%is given by this number.
%-Rewards: Array which contains the
%rewards calculated for each solution returned. Since we are using log terms in
%the assignment matrices, taking the exponential of these would give
%probability factors. 


[NofCustomers,NofItems]=size(RewardMatrix);%Get the number of customers and items

Customer2Item=zeros(N,NofCustomers);%Define the output array for customer assignments
Item2Customer=zeros(N,NofItems);%Define the output array for item assignments (Assignment 0 means no assignment)
Rewards=zeros(N,1);%Define the reward values for the solutions

[FirstCustomer2Item,FirstItem2Customer]=Auction(RewardMatrix);%Get the best solutions


%Murty Algorithm Problem-Solution pairs list
DynamicRewardMatrix{1}=RewardMatrix;%Assignment matrix of the Problems in the List
DynamicCustomer2Item=FirstCustomer2Item;%Customer Assignments for the Solutions in the List
DynamicItem2Customer=FirstItem2Customer;%Item Assignments for the Solutions in the List
DynamicReward=RewardCalculate(RewardMatrix,FirstCustomer2Item);%Rewards for the Solutions in the List  


minval=min(RewardMatrix(:))-10^7;%Instrumental Minimum value to Handle Solution Censoring

i=0;%Set the index to zero.
while  (i<=N) && ~isempty(DynamicReward),
    i=i+1; %increment i. 
    %Find the maximum reward solution in the list
    [maxval indmax]=max(DynamicReward);%Find the maximum reward solution in the list   
    TempCustomer2Item=DynamicCustomer2Item(indmax,:);%get the customer assignments of the solution
    TempItem2Customer=DynamicItem2Customer(indmax,:);%get the item assignments of the solution
    TempRewardMatrix=DynamicRewardMatrix{indmax};%get the rewardmatrix of the problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Add the maximum reward solution in the list to the list of solutions
    %to be returned
    Customer2Item(i,:)=TempCustomer2Item;%Set the ith best solution's customer assignments
    Item2Customer(i,:)=TempItem2Customer;%Set the ith best solution's item assignments
    Rewards(i)=maxval;%Set the ith best solution's reward
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Remove the maximum reward solution from the list
    DynamicReward(indmax)=[];%Remove the reward from the list
    DynamicCustomer2Item(indmax,:)=[];%Remove the customer assignments from the list
    DynamicItem2Customer(indmax,:)=[];%Remove the item assignments from thel list
    DynamicRewardMatrix(indmax)=[];%Remove the assignment matrix of the problem from the list
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (i==N)%if we have already determined the Nth best solution
        break;%no need to continue to spend computation time, get out of the loop.
    end
    for j=1:NofCustomers,%Trace all associations of the given solution
        TempRewardMatrix(j,TempCustomer2Item(j))=minval;%Make the customer association invalid
        [LoopCustomer2Item,LoopItem2Customer]=Auction(TempRewardMatrix);%solve the modified problem
        ValidInvalidFlag=CheckValidity(TempRewardMatrix,LoopCustomer2Item,minval);%check the validity of the new solution            
        if ValidInvalidFlag %if the found solution is valid
           Ndynamic=length(DynamicReward);%find the number of problem-solution pains in the list
           DynamicReward(Ndynamic+1)=RewardCalculate(TempRewardMatrix,LoopCustomer2Item);%Add the reward of thesolution to the list
           DynamicCustomer2Item(Ndynamic+1,:)=LoopCustomer2Item;%Add customer assignments to the list
           DynamicItem2Customer(Ndynamic+1,:)=LoopItem2Customer;%Add the item assignments to the list 
           DynamicRewardMatrix{Ndynamic+1}=TempRewardMatrix;%Add the assignment matrix of the problem to the list
        end
        %restore the reward value
        TempRewardMatrix(j,TempCustomer2Item(j))=RewardMatrix(j,TempCustomer2Item(j));%Make the customer association valid again
        %remove all possible associations except the selected one
        %i.e., make all association except the customer association invalid
        TempRewardMatrix(j,1:(TempCustomer2Item(j)-1))=minval;
        TempRewardMatrix(j,(TempCustomer2Item(j)+1):end)=minval;
        TempRewardMatrix(1:(j-1),TempCustomer2Item(j))=minval;
        TempRewardMatrix((j+1):end,TempCustomer2Item(j))=minval;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

Nsolutions=i;%Find the number of solutions obtained
if i~=N,%if we did not get required number of solutions
   %Discard the unused elements of the output arrays 
   Customer2Item=Customer2Item(1:i,:) ;
   Item2Customer=Item2Customer(1:i,:) ;
   Rewards=Rewards(1:i);
end


function Reward=RewardCalculate(RewardMatrix,Customer2Item)
%Calculate the reward of a Customer Assignment  
mysum=0;%Set reward to zero
NofCustomers=length(Customer2Item);
for i=1:NofCustomers%For all customers
    mysum=mysum+RewardMatrix(i,Customer2Item(i));%sum the reward of thwe item assigned to the customer
end
Reward=mysum;%return the sum

function ValidInvalidFlag=CheckValidity(RewardMatrix,Customer2Item,minval)
%Calculate the validity of a customer assignment
%A customer assignment to an item with reward of minval 
%represents an infeasible assignment. This means that 
%such solutions are invalid
ValidInvalidFlag=1;%valid at the beginning
NofCustomers=length(Customer2Item);
for i=1:NofCustomers
    if (RewardMatrix(i,Customer2Item(i))==minval)%if we have minval reward association
        ValidInvalidFlag=0;%make the solution invalid
    end
end

