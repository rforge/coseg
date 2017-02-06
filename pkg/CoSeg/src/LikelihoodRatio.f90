!This fortran code calculates the likelihood ratio presented in the method in the Mohammadi paper "A simple method for cosegregation analysis to evaluate the pathogenicity of unclassified variants; BRCA1 and BRCA2 as an example"



! module constants
!
! 	implicit none
!
! 	integer, parameter :: dble_prec = kind(0.0d0)
! 	integer, parameter :: real_prec = kind(0.0e0)
! 	real(kind=dble_prec), parameter :: pi = 3.14159
! 	real(kind=dble_prec), parameter :: zero = 0.00000000
! 	real(kind=dble_prec), parameter :: one = 1.00000000
! 	real(kind=dble_prec), parameter :: two = 2.00000000
! 	real(kind=dble_prec), parameter :: ten = 10.0000000
! 	!	real(kind=dble_prec), parameter :: root2=sqrt(2)
! 	real(kind=dble_prec), parameter :: eps=1.e-10
! 	real(kind=dble_prec), parameter :: delta=one/16.000!1.e-2
! 	real(kind=dble_prec), parameter :: pseudo_count=1.0000/10.0000
!
!
! end module constants



! module rprint
!
! 	use constants
!
! 	implicit none
!
! 	contains
!
!
! ! 	subroutine dblepr(line,nchar,dble_vector,ndata)
! ! 	!line is a char up to 255
! ! 	!nchar is line's length or -1 if whole line is to be printed
! ! 	!dble_vector is vector to be printed of length at least ndata
! ! 	!ndata is the number of elements of vector to print
!
! ! 	implicit none
!
! ! 	character(len=*) :: line
! ! 	integer :: nchar, ndata
! ! !	real(kind=dble_prec), dimension(:) ::dble_vector
! ! 	real(kind=dble_prec) :: dble_vector
!
! ! 	if(nchar==-1) then
! ! 		nchar=len(trim(line))
! ! 	end if
!
! ! 	print*, line(1:nchar)
! ! !	print*,dble_vector(1:ndata)
! ! 	print*,dble_vector
!
! ! 	end subroutine dblepr
!
!
!
!
! ! 	subroutine intpr(line,nchar,int_vector,ndata)
! ! 	!line is a char up to 255
! ! 	!nchar is line's length or -1 if whole line is to be printed
! ! 	!int_vector is vector to be printed of length at least ndata
! ! 	!ndata is the number of elements of vector to print
!
! ! 	implicit none
!
! ! 	character(len=*) :: line
! ! 	integer :: nchar, ndata
! ! 	integer ::int_vector
!
! ! 	if(nchar==-1) then
! ! 		nchar=len(trim(line))
! ! 	end if
!
! ! 	print*, line(1:nchar)
! ! !	print*,int_vector(1:ndata)
! ! 	print*,int_vector
!
! ! 	end subroutine intpr
!
!
!
!
!
! end module rprint



!
!
! subroutine calculate_likelihood_ratio_fast(pedid, pedmomid, peddadid, pedfemale, &
!   PedGenotype, pedage, pedaffected, NumberPeople, LikelihoodRatio, separating_meioses)
!
! use rprint
! use constants
!
! implicit none
!
! !here we list the needed inputs.
! integer :: NumberPeople
! integer, dimension(NumberPeople) :: pedid, pedmomid, peddadid,  PedGenotype, pedaffected
! real(kind=dble_prec), dimension(NumberPeople) :: pedage
! logical, dimension(NumberPeople) :: pedfemale
!
! !outputs since R does not like functions
! real(kind=dble_prec) :: LikelihoodRatio
! integer :: separating_meioses
!
! !here are the other variables
! real(kind=dble_prec) :: lr_numerator, lr_denominator, check_num_genotype_probability, check_den_genotype_probability
! integer :: i,j,k,number_genotypes_found, counter, int1
! integer, dimension(NumberPeople) :: momrow, dadrow, degree
! logical :: logical1, logical2
! logical, dimension(NumberPeople) :: observed_genotype, logical_vec, pedigree_founder, &
! 						ancestor_vec, future_vec, current_vec
! logical, dimension(NumberPeople,NumberPeople) :: AncestorDescendentArray
!
!
!
!
! !adding vectors that contain the row numbers of the parents instead of ids
! momrow=0
! dadrow=0
! do i=1,NumberPeople
!   if(pedmomid(i).ne.0) then
!     logical1=.false.
!     logical2=.false.
!     do j=1,NumberPeople
!       if(.not. logical1) then
!         if(pedmomid(i).eq.pedid(j)) then
!           momrow(i)=j
!           logical1=.true.
!         end if
!       end if
!       if(.not. logical2) then
!         if(peddadid(i).eq.pedid(j)) then
!           dadrow(i)=j
!           logical12=.true.
!         end if
!       end if
!       if(logical1 .and. logical2) then
!         exit
!       end if
!     end do
!   end if
! end do
!
! !adding a boolean vector that shows which individuals have an observed genotype
! observed_genotype=.false.
! do i=1,NumberPeople
!   if(PedGenotype(i).ne.2) then
!     observed_genotype(i)=.true.
!   end if
! end do
!
! !adding in degree information
! !start with first row individual.  Assign degree 1.  Assign all offspring degree 2 and parents degree 0 if any.  Move to all newly assigned individuals and assign degrees that are either +-1 of the current degree of the individual based off parent/offspring relationship.  In the end, add a constant to the degree vector to make the lowest degree 1.
! degree=0
! degree(1)=1
! logical_vec=.false.
! do while(.not. all(logical_vec))
!   do i=1,NumberPeople
!     if(.not. logical_vec(i) .and. degree(i).ne.0) then
!       !sets parents to degree -1
!       if(momrow(i).ne.0) then
!         degree(momrow(i))=degree(i)-1
!         degree(dadrow(i))=degree(i)-1
!       end if
!       !sets offspring to degree +1
!       if(pedfemale(i)) then
!         do j=1,NumberPeople
!           if(momrow(j).eq.i) then
!             degree(j)=degree(i)+1
!           end if
!         end do
!       else
!         do j=1,NumberPeople
!           if(dadrow(j).eq.i) then
!             degree(j)=degree(i)+1
!           end if
!         end do
!       end if
!
!       logical_vec(i)=.true.
!     end if
!   end do
! end do
! degree=degree+abs(minval(degree))+1
!
!
! !Here, we find all ancestors and descendents of each individual
! !Note that ancestor.descendent.array[i,j]=TRUE if j is an ancestor of i or i is a descendent of j
! !Also note that for convenience in the code, ancestor.descendent.array[i,i]=TRUE
!
! AncestorDescendentArray=.false.
! pedigree_founder=.false.
! do i=1,NumberPeople
!   if(momrow(i).eq.0) then
!     pedigree_founder(i)=.true.
!   end if
! end do
! do i=1,NumberPeople
!   ancestor_vec=.false.
!   future_vec=.false.
!   ancestor_vec(i)=.true.
!   current_vec=ancestor_vec
!   logical1=.true.
!   do while(logical1)
!     do j=1,NumberPeople
!       if(current_vec(j) .and. .not. pedigree_founder(j)) then
!         future_vec(dadrow(j))=.true.
!         future_vec(momrow(j))=.true.
!       end if
!     end do
!     if(any(futurevec)) then
!       ancestor_vec=ancestor_vec .or. future_vec
!       current_vec=future_vec
!       future_vec=.false.
!     else
!       logical1=.false.
!       exit
!     end if
!   end do
!   AncestorDescendentArray(i,:)=ancestor_vec
! end do
!
!
! !The methods below depend on the ordering of the pedigree.  The ancestors of the pedigree should have a lower number than the descendents of the pedigree.  Here, we check to see if that's the case and reorder if not.
! counter=1
! int1=1
!
!
! !change of plans... lets start with recoding the crux of the function since that is what runs for >99.99% of the time...
!
!
! end subroutine calculate_likelihood_ratio_fast
!
!





subroutine likelihood_ratio_main(NumberPeople, NumberProbandFounders, &
	NumberPossibleFounders, ObservedSeparatingMeioses, NumberOffspring, &
	PedGenotype, FounderCols, &
	NumberGenotypesVec, AllPhenotypeProbabilities, ProbandAncestors,  &
	ObservedVector, AncestorDescendentArray, MomRow, DadRow, &
	MinimalObservedPedigree, LikelihoodRatio)
  !the general idea here is that we are finding all possible genotype.  The way we do this is we start with the left most variable genotype(lowest number) and fix it to 0.  This in turn fixes a lot of the genotype to the right(higher numbers).  To keep track of these, we set status.vector to be an integer vector where 0 means the genotype does not need to be modified anymore, otherwise it shows the number of the genotype that currently fixed it or number.people+1 if it hasn't been touched yet.  We then go right(increasing number) to the next variable genotype and fix it to 0, fix those that become fixed because of that, and modify status.vector accordingly.  Once status.vector is completely fixed, we have a viable genotype.  We now save the genotype*phenotype probabilities for that genotype for the numerator and denominator of the likelihood ratio.  We then proceed back left left(decreasing number) to the last fixed genotype and fix it to 1 and see if this fixes anything to the right(though it shouldn't).

! use constants
! use rprint

implicit none

!constants
integer, parameter :: dble_prec = kind(0.0d0)

!inputs
integer :: NumberPeople, NumberProbandFounders, NumberPossibleFounders, &
	ObservedSeparatingMeioses
integer, dimension(NumberPeople) :: NumberOffspring, PedGenotype, MomRow, DadRow
integer, dimension(NumberProbandFounders) :: FounderCols
integer, dimension(2) :: NumberGenotypesVec
real(kind=dble_prec), dimension(2,NumberPeople) :: AllPhenotypeProbabilities
logical, dimension(NumberPeople) :: MinimalObservedPedigree, ProbandAncestors, &
	ObservedVector
logical, dimension(NumberPeople,NumberPeople) :: AncestorDescendentArray

!outputs
real(kind=dble_prec) :: LikelihoodRatio

!extra variables
integer :: i, j, k, counter, int1, int2
integer :: number_genotypes_found
integer, dimension(NumberPeople) :: status_vector, temp_changes
logical, dimension(NumberPeople) :: founder_descendents, lineage_genotype, &
  variable_genotype, temp_genotype, logical_vec!, implied_carriers
logical :: logical1
real(kind=dble_prec) :: phenotype_probability, check_num_genotype_probability, &
  check_den_genotype_probability, denominator_genotype_probability, &
  numerator_genotype_probability, lr_denominator, lr_numerator, tolerance
real(kind=dble_prec), dimension(NumberPeople) :: phenotype_probabilities_ratios, &
	negative_powers_of_two



!initializing some variables
NumberGenotypesVec=0
phenotype_probabilities_ratios=AllPhenotypeProbabilities(1,:)/AllPhenotypeProbabilities(2,:)
number_genotypes_found=0
lr_denominator=0
lr_numerator=0
check_den_genotype_probability=0
check_num_genotype_probability=0

!storing all useful negative powers of two for analysis
do i=1,NumberPeople
	negative_powers_of_two(i)=2.**(-i)
end do

do i=1,NumberProbandFounders
  !at least 2 of the proband.founders are interchangeable so we don't need to enumerate through them separately.  We could improve performance drastically by finding any interchangeable proband.founders and intelligently cycling through them.

  !link up the founder and the proband.  Figure out which genotype must be fixed between the current founder and the proband
  !starts at the current proband founder and finds all possible descendents
  founder_descendents=AncestorDescendentArray(:,FounderCols(i))

  !Assumes that the intersection of all ancestors of the proband with all descendents of the current founder give the lineage from the founder to the proband.  Note that these must be fixed to 1 or carrier
  lineage_genotype=ProbandAncestors .and. founder_descendents
  !these genotype are the only ones that can be 0 or 1
  variable_genotype=founder_descendents .and. .not.lineage_genotype

	!A vector of the union of the intersection of the current founder's descendents and the all known carriers ancestors
	! implied_carriers=.FALSE.
	! do j=1,NumberPeople
	! 	if(PedGenotype(j).eq.1) then
	! 		implied_carriers=implied_carriers .or. (founder_descendents .and. AncestorDescendentArray(j,:))
	! 	end if
	! end do
  !Here we go need to go through variable.genotype and starting from the one end set it to 0 or 1 and fix all the genotype affected by that choice, fix the next unfixed one, and so on until all are fixed.  We then work backwards to get all the possible choices.
  !Since I'm starting from 1, assume descendents outnumber ancestors meaning that being a non-carrier would fix more genotype than being a carrier.  Note this assumption shouldn't affect the final result.  In other words, start assuming j is a non-carrier... temp.genotype[j]=0
  temp_genotype=lineage_genotype !temp.genotype is a boolean vector that holds the current individual's carrier status.

  do j=1,NumberPeople
    if(variable_genotype(j)) then
      status_vector(j)=NumberPeople+1
    else
      status_vector(j)=0
    end if
  end do

  phenotype_probability=1
  do j=1,NumberPeople
    if(temp_genotype(j)) then
      phenotype_probability=phenotype_probability*AllPhenotypeProbabilities(2,j)
    else
      phenotype_probability=phenotype_probability*AllPhenotypeProbabilities(1,j)
    end if
  end do

  do while(sum(status_vector)>=0)

    !note: temp.genotype is 0 everywhere except for lineage.genotype
    !move right(increasing) filling the vector
    do j=1,NumberPeople
      if(status_vector(j).eq.NumberPeople+1) then
        !Assume j is set to 0.  Set all descendents of j to be fixed at 0
        do k=1,NumberPeople
          if(AncestorDescendentArray(k,j)) then
            status_vector(k)=j
          end if
        end do
      end if
    end do

    !found a viable vector so save genotype*phenotype probabilities into numerator and denominator
		if(number_genotypes_found==2147483647) then
			number_genotypes_found=0
			NumberGenotypesVec(1)=NumberGenotypesVec(1)+1
		else
    	number_genotypes_found=number_genotypes_found+1
		end if
		!call intpr("temp_genotype",-1,temp_genotype,NumberPeople)

		!prints out the current status
		if(mod(number_genotypes_found,25000000).eq.0) then
			call intpr("Number genotypes found times 25 million: ",-1,number_genotypes_found/25000000,1)
			call dblepr("Fraction done (note that this is note linear with number genotypes):",-1,check_den_genotype_probability,1)
		end if

    !denominator is updated regardless
		!int1 is the number of offspring that are potential carriers.
    int1=0
		! if(NumberProbandFounders.eq.1) then
		! 	int1=-1
		! end if
    do j=1,NumberPeople
      if(temp_genotype(j)) then
        int1=int1+NumberOffspring(j)
      end if
    end do

		if(int1>0 .and. int1<=NumberPeople) then
			denominator_genotype_probability=negative_powers_of_two(int1)
		else
    	denominator_genotype_probability=2.**(-int1)
		end if
    lr_denominator=lr_denominator+phenotype_probability*denominator_genotype_probability
    check_den_genotype_probability=check_den_genotype_probability+denominator_genotype_probability

    !numerator is updated if the current genotype vector matches the observed genotype
    logical1=.true.
    do j=1,NumberPeople
      if(ObservedVector(j)) then
        if(temp_genotype(j)) then
          if(PedGenotype(j).eq.0) then
            logical1=.false.
            exit
          end if
        else !if temp_genotype(j).eq.0
          if(PedGenotype(j).eq.1) then
            logical1=.false.
            exit
          end if
        end if
      end if
    end do
    if(logical1) then
			if(all(ObservedVector)) then
				numerator_genotype_probability=1
			else

				!Here we count the number of unknown genotype people who have a carrier parent
				!note that implied genotypes are counted as known (Thus MinimalObservedPedigree is used)
				int2=1
				if(NumberPossibleFounders.eq.1) then
					int2=0
				end if
		    do j=1,NumberPeople
					if(.not.ObservedVector(j) .and. .not.MinimalObservedPedigree(j)) then !check his parents if one of them is a carrier
						if(MomRow(j)>0 .and. DadRow(j)>0) then !check if the person is a founder
							if(temp_genotype(MomRow(j)) .or. temp_genotype(DadRow(j))) then
			        	int2=int2+1
							end if
						end if
		      end if
		    end do

				if(int2>0 .and. int2<=NumberPeople) then
					numerator_genotype_probability=negative_powers_of_two(int2)
				! call dblepr("1st Numerator check: ",-1,numerator_genotype_probability,1)
			!else if(ObservedSeparatingMeioses-int1>0) then
      	!numerator_genotype_probability=2.**(-ObservedSeparatingMeioses+int1)
				else
					numerator_genotype_probability=2.**(-int2)
					! call dblepr("Numerator check: ",-1,numerator_genotype_probability,1)
				end if
			end if
      check_num_genotype_probability=check_num_genotype_probability+numerator_genotype_probability
      lr_numerator=lr_numerator+phenotype_probability*numerator_genotype_probability
    end if

    !move left(decreasing) and change last move to have genotype 1.
    int1=maxval(status_vector)
    if(int1>0) then
      do j=1,NumberPeople
        if(status_vector(j).eq.int1) then
          logical_vec(j)=.true.
        else
          logical_vec(j)=.false.
        end if
      end do
      temp_changes=0
      counter=0
      do j=1,NumberPeople
        if(logical_vec(j).and.temp_genotype(j)) then
          counter=counter+1
          temp_changes(counter)=j
        end if
      end do
      if(counter>0) then
        do j=1,counter
          phenotype_probability=phenotype_probability*phenotype_probabilities_ratios(temp_changes(j))
        end do
      end if
      do j=1,NumberPeople
        if(logical_vec(j)) then
          temp_genotype(j)=.false.
          status_vector(j)=NumberPeople+1
        end if
      end do
      phenotype_probability=phenotype_probability/phenotype_probabilities_ratios(int1)
      temp_genotype(int1)=.true.
      status_vector(int1)=maxval(status_vector(1:(int1-1)))
    else
      status_vector(1)=-1
    end if
  end do
end do

LikelihoodRatio=lr_numerator/lr_denominator

NumberGenotypesVec(2)=number_genotypes_found

tolerance=1E-6

if(check_num_genotype_probability<1-tolerance .OR. check_num_genotype_probability>1+tolerance) then
  call dblepr("Error: Total numerator genotype probability(Should be 1): ",-1,check_num_genotype_probability,1)
end if
if(check_den_genotype_probability<1-tolerance .OR. check_den_genotype_probability>1+tolerance) then
  call dblepr("Error: Total denominator genotype probability(Should be 1): ",-1, check_den_genotype_probability,1)
end if

!call intpr("Number of genotypes found: ",-1,number_genotypes_found,1)
!call dblepr("Likelihood ratio: ",-1,LikelihoodRatio,1)


end subroutine likelihood_ratio_main
















!end of file
