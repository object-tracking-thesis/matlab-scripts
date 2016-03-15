function [Customer2Item,Item2Customer]=Auction(RewardMatrix)



%Auction algorithm
%based on the description in the book by Blackman & Popoli (1999) 
%Solves the assignment problem (maximization)

%Input: RewardMatrix
%Requires number of columns to be greater than or equal to the number of
%rows

%Output: Customer2Items, Item2Customer
%-The auction algorithm is in analogy with an actual auction.
%-The rows of RewardMatrix corresponds to customers
%-The columns of RewardMatrix corresponds to items in the auction
%-This auction algorithm requires the number of items greater than or equal
%to the customers. In tracking, the assignment matrices (Slide 18 of Lecture 4)  
%we form always satisfy this.
%-The algorithm assigns each customer an item.
%-Some items might not be assigned to any customer.

%Customer2Item is a vector of size equal to the number of customers (measurements in our case)
%Customer2Item(i) is the index of the item (track) assigned to the
%ith customer (measurement)

%Item2Customer is a vector of size equal to the number of columns of
%RewardMatrix (generally this is sum of the number of targets + number of measurements). 
%Item2Customer(i) is the index of the customer which assigned to the
%ith item. If Item2Customer(i)==0, this means the item is not assigned.

epsilon=0.1; %amount of deviation from the optimal reward

[NofCustomers,NofItems]=size(RewardMatrix);
if (NofCustomers>NofItems)
    error('Number of columns must be greater than or equal to the number of rows');
end

Item2Customer=zeros(1,NofItems);
Customer2Item=zeros(1,NofCustomers);

while ~isempty(find(Customer2Item==0)),
    if (NofItems==1) %if there is only one item
        [maxval Item2Customer]=max(RewardMatrix);%Assign the item to the best customer
        Customer2Item(Item2Customer)=1;%Assign the corresponding customer to the item
    else
        for i=1:NofCustomers,
            if ~Customer2Item(i),
                [maxval,maxind]=max(RewardMatrix(i,:));%find maximum element value and its index
                RewardMatrix(i,maxind)=min(RewardMatrix(i,:))-1;%make the maximum minimum to find second maximum
                [secondmaxval,secondmaxind]=max(RewardMatrix(i,:));%find the second maximum value and its index
                RewardMatrix(i,maxind)=maxval; %restore the maximum value

                Customer2Item(i)=maxind; %Assign the customer the item
                if Item2Customer(maxind),%if item is already assigned
                    Customer2Item(Item2Customer(maxind))=0;%unassign the corresponding customer
                end
                Item2Customer(maxind)=i; %Assign the item to the customer
                RewardMatrix(:,maxind)=RewardMatrix(:,maxind)-(maxval-secondmaxval+epsilon);%reduce the item's value
            end
        end
    end
end


